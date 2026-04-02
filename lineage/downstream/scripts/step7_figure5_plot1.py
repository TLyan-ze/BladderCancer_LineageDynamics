import os
import sys


from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.spatial import distance
from scipy.cluster import hierarchy
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from tqdm.auto import tqdm
sys.path.append("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/KPTracer-release-main/cassiopeia-kp/") # specify path to jungle
import cassiopeia
from cassiopeia.Analysis import small_parsimony
sys.path.append("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/plot/Figure5_S5/scripts/")
import tree_utilities


# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))

# 将当前工作目录设置为脚本所在目录
os.chdir(script_dir)

##########################################################################################################

adata = sc.read_h5ad(f'/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/expression/scRNA4_3000_25.h5ad')

tumor_list = pd.read_csv(f"{script_dir}/tumor_list.txt", sep='\t')['Tumor'].values

tumors = [
    t
    for t in tumor_list
    if "Fam" not in t
    #and "Met" not in t
    and "All" not in t
    #and 'NT' in t
    #and t.split("_")[2].startswith("T")
]
print(tumors)
corr_dict = {}

for tumor in tqdm(tumors):
    print (tumor)
##############################################################################################################################
    #tumor = 'Met_MGH_2'

    graph = tree_utilities.prepare_tumor_tree(tumor, adata,tree_dir = f"/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele_new/trees/",
                                                FILTER_PROP=0.025,
                                                column='leiden_sub')

    for u, v in graph.edges():
        _length = graph[u][v]["length"]
        if _length > 0:
            _length = 1
        graph[u][v]["length"] = _length

    phylogenetic_distance_matrix, leaf_pairs, edit_distance_matrix, tree_diameter, n_targets = tree_utilities.compute_pairwise_dist_nx(graph)

    _leaves = [n for n in graph if graph.out_degree(n) == 0]
    np.fill_diagonal(phylogenetic_distance_matrix.values, 0)
    leaf_states = adata.obs.loc[_leaves, 'leiden_sub']

    observed_inter_cluster_df = tree_utilities.get_inter_cluster_df(leaf_states, phylogenetic_distance_matrix)

    B = 1000
    background = defaultdict(list)
    for _ in tqdm(range(B)):
        permuted_assignments = leaf_states.copy()
        permuted_assignments.index = np.random.permutation(leaf_states.index.values)
        bg_df = tree_utilities.get_inter_cluster_df(permuted_assignments, phylogenetic_distance_matrix)
        
        for s1 in bg_df.index:
            for s2 in bg_df.index:
                background[(s1, s2)].append(bg_df.loc[s1, s2])
                
    null_means = observed_inter_cluster_df.copy()
    null_sds = observed_inter_cluster_df.copy()

    for s1 in null_means.index:
        for s2 in null_means.columns:
            null_means.loc[s1, s2] = np.mean(background[(s1, s2)])
            null_sds.loc[s1, s2] = np.std(background[(s1, s2)])


    z_df = observed_inter_cluster_df.copy()
    for ind in z_df.index:
        for col in z_df.columns:
            z_df.loc[ind, col] = (z_df.loc[ind, col] - null_means.loc[ind, col]) / null_sds.loc[ind, col]
    z_df.fillna(0, inplace=True)

    mu = np.mean(z_df.values.ravel())
    sigma = np.std(z_df.values.ravel())

    zz_df = z_df.apply(lambda x: (x - mu) / sigma, axis=1)



    uniq_labels = adata.obs['leiden_sub'].astype(int).unique()
    _normalized_cluster_df = pd.DataFrame(np.zeros((len(uniq_labels), len(uniq_labels))), index=uniq_labels, columns = uniq_labels)
    for i in zz_df.index:
        for j in zz_df.columns:
            _normalized_cluster_df.loc[i, j] = (zz_df.max().max() - zz_df).loc[i, j]

    tree_utilities.plot_graph_on_umap(adata, _normalized_cluster_df, weight=1,
                                    cluster_column = 'leiden_sub', title=f'Evolutionary Coupling, {tumor}')
    plt.savefig("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele_new/Figure/Figure5_S5/plot/Evolutionary_Coupling"+str(tumor)+".png")
    plt.savefig("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele_new/Figure/Figure5_S5/plot/Evolutionary_Coupling"+str(tumor)+".pdf")

