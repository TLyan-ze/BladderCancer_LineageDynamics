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




##########################################################################################################
data_directory = "/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/"

adata = sc.read_h5ad(f'{data_directory}/expression/scRNA4.h5ad')

tumor_list = pd.read_csv(f"{data_directory}/plot/Figure5_S5/tumor_list.txt", sep='\t')['Tumor'].values

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

    graph = tree_utilities.prepare_tumor_tree(tumor, adata,tree_dir = f"/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele/trees/",
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

    df=_normalized_cluster_df
    min_prop = 0.0
    weight=1.0
    cluster_column='leiden_sub'
    title=f'Evolutionary Coupling, {tumor}'
    umap_coords = pd.DataFrame(adata.obsm['X_umap'], index = adata.obs_names)
    paga_pos = {}
    for n, g in adata.obs.groupby(cluster_column):

        paga_pos[int(n)] = umap_coords.loc[g.index].mean(0).values

    adjacency_solid = df.fillna(0).copy()
    adjacency_solid.index = [int(n) for n in adjacency_solid.index]
    adjacency_solid.columns = [int(n) for n in adjacency_solid.columns]
    adjacency_solid = adjacency_solid.loc[adjacency_solid.index, adjacency_solid.columns]
    import matplotlib as mpl
    import sys
    import ete3
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import networkx as nx
    import numba
    import numpy as np
    import pandas as pd
    import scanpy as sc
    adjacency_solid[adjacency_solid < min_prop] = 0.0

    base_edge_width = weight * mpl.rcParams['lines.linewidth']

    nx_g_solid = nx.Graph(adjacency_solid.values)

    widths = [x[-1]['weight'] for x in nx_g_solid.edges(data=True)]
    widths = base_edge_width * np.array(widths)

    h = plt.figure(figsize=(7,7))
    ax = plt.gca()
    nx.draw_networkx_nodes(nx_g_solid, paga_pos, ax=ax)
    nx.draw_networkx_labels(nx_g_solid, paga_pos, ax=ax)
    nx.draw_networkx_edges(
                    nx_g_solid, paga_pos, ax=ax, width=widths, edge_color='black'
                )
    sc.pl.embedding(adata, basis='X_umap', color=cluster_column, ax=ax, show=False)
    if title:
        plt.title(title)
    plt.tight_layout()
    ax.set_xticks([])
    ax.set_yticks([])
    plt.show()
    #图一
    plt.savefig("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/plot/Figure5_S5/plot/Evolutionary_Coupling"+str(tumor)+".png")
    plt.savefig("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/plot/Figure5_S5/plot/Evolutionary_Coupling"+str(tumor)+".pdf")

################################################################################################################################

    K = 10

    graph = tree_utilities.prepare_tumor_tree(tumor, adata,
                                            tree_dir = f"{data_directory}/allele/trees",
                                            FILTER_PROP=0.025,
                                            column='leiden_sub')

    phylogenetic_distance_matrix = pd.read_csv(f'./data/{tumor}_phylogenetic_distance_matrix.tsv', sep='\t', index_col = 0)
    edit_distance_matrix = pd.read_csv(f'./data/{tumor}_edit_distance_matrix.tsv', sep='\t', index_col = 0)

    _leaves = np.intersect1d(phylogenetic_distance_matrix.index.values, [n for n in graph if graph.out_degree(n) == 0])
    _leaves = np.intersect1d(_leaves, adata.obs_names)

    leaf_states = adata.obs.loc[_leaves, 'leiden_sub']

    phylogenetic_distance_matrix = phylogenetic_distance_matrix.loc[_leaves, _leaves]
    edit_distance_matrix = edit_distance_matrix.loc[_leaves, _leaves]

    phylo_inter_cluster_df = tree_utilities.get_inter_cluster_df(leaf_states, phylogenetic_distance_matrix)
    edit_inter_cluster_df = tree_utilities.get_inter_cluster_df(leaf_states, edit_distance_matrix)
    kth_inter_cluster_df = tree_utilities.get_inter_cluster_df(leaf_states, phylogenetic_distance_matrix, tree_utilities.average_nn_dist, k = K)

    B = 500
    phylo_background = defaultdict(list)
    edit_background = defaultdict(list)
    kth_background = defaultdict(list)

    for _ in tqdm(range(B)):
        permuted_assignments = leaf_states.copy()
        permuted_assignments.index = np.random.permutation(leaf_states.index.values)
        bg_df = tree_utilities.get_inter_cluster_df(permuted_assignments, phylogenetic_distance_matrix)
        
        for s1 in bg_df.index:
            for s2 in bg_df.index:
                phylo_background[(s1, s2)].append(bg_df.loc[s1, s2])
                
        bg_df = tree_utilities.get_inter_cluster_df(permuted_assignments, edit_distance_matrix)
        
        for s1 in bg_df.index:
            for s2 in bg_df.index:
                edit_background[(s1, s2)].append(bg_df.loc[s1, s2])
                
        bg_df = tree_utilities.get_inter_cluster_df(permuted_assignments, phylogenetic_distance_matrix, func=tree_utilities.average_nn_dist, k=K)
        
        for s1 in bg_df.index:
            for s2 in bg_df.index:
                kth_background[(s1, s2)].append(bg_df.loc[s1, s2])
                
    phylo_null_means = phylo_inter_cluster_df.copy()
    phylo_null_sds = phylo_inter_cluster_df.copy()

    edit_null_means = edit_inter_cluster_df.copy()
    edit_null_sds = edit_inter_cluster_df.copy()

    kth_null_means = kth_inter_cluster_df.copy()
    kth_null_sds = kth_inter_cluster_df.copy()


    for s1 in phylo_null_means.index:
        for s2 in phylo_null_means.columns:
            phylo_null_means.loc[s1, s2] = np.mean(phylo_background[(s1, s2)])
            phylo_null_sds.loc[s1, s2] = np.std(phylo_background[(s1, s2)])
            
            edit_null_means.loc[s1, s2] = np.mean(edit_background[(s1, s2)])
            edit_null_sds.loc[s1, s2] = np.std(edit_background[(s1, s2)])
            
            kth_null_means.loc[s1, s2] = np.mean(kth_background[(s1, s2)])
            kth_null_sds.loc[s1, s2] = np.std(kth_background[(s1, s2)])
    phylo_zscores = phylo_inter_cluster_df.copy()
    edit_zscores = edit_inter_cluster_df.copy()
    kth_zscores = kth_inter_cluster_df.copy()

    for ind in phylo_zscores.index:
        for col in phylo_zscores.columns:
            phylo_zscores.loc[ind, col] = (phylo_zscores.loc[ind, col] - phylo_null_means.loc[ind, col]) / phylo_null_sds.loc[ind, col]
            edit_zscores.loc[ind, col] = (edit_zscores.loc[ind, col] - edit_null_means.loc[ind, col]) / edit_null_sds.loc[ind, col]
            kth_zscores.loc[ind, col] = (kth_zscores.loc[ind, col] - kth_null_means.loc[ind, col]) / kth_null_sds.loc[ind, col]
            
    phylo_zscores.fillna(0, inplace=True)
    edit_zscores.fillna(0, inplace=True)
    kth_zscores.fillna(0, inplace=True)

    mu = np.mean(phylo_zscores.values.ravel())
    sigma = np.std(phylo_zscores.values.ravel())
    phylo_zscores_norm = phylo_zscores.apply(lambda x: (x - mu) / sigma, axis=1)

    mu = np.mean(edit_zscores.values.ravel())
    sigma = np.std(edit_zscores.values.ravel())
    edit_zscores_norm = edit_zscores.apply(lambda x: (x - mu) / sigma, axis=1)

    mu = np.mean(kth_zscores.values.ravel())
    sigma = np.std(kth_zscores.values.ravel())
    kth_zscores_norm = kth_zscores.apply(lambda x: (x - mu) / sigma, axis=1)


    pscores = np.exp(-phylo_zscores_norm / phylo_zscores_norm.abs().max().max())
    escores = np.exp(-edit_zscores_norm / edit_zscores_norm.abs().max().max())
    kscores = np.exp(-kth_zscores_norm / kth_zscores_norm.abs().max().max())

    corr_edit = round(stats.pearsonr(pscores.values.ravel(), escores.values.ravel())[0],4)

    h = plt.figure(figsize=(5,5))
    plt.scatter(pscores.loc[edit_zscores_norm.index, edit_zscores_norm.columns].values.ravel(), escores.values.ravel(), color='black')
    plt.plot(pscores.values.ravel(), pscores.values.ravel(), color='red')
    plt.xlabel("Phylogenetic EvoCoupling")
    plt.ylabel("Edit EvoCoupling")
    plt.title(f"{tumor}, corr: {corr_edit}")

    plt.savefig("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/plot/Figure5_S5/plot/edit"+str(tumor)+".png")
    plt.savefig("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/plot/Figure5_S5/plot/edit"+str(tumor)+".pdf")

    corr_kth = round(stats.pearsonr(pscores.values.ravel(), kscores.values.ravel())[0],4)

    h = plt.figure(figsize=(5,5))
    plt.scatter(pscores.loc[kscores.index, kscores.columns].values.ravel(), kscores.values.ravel(), color='black')
    plt.plot(pscores.values.ravel(), pscores.values.ravel(), color='red')
    plt.xlabel("Phylogenetic EvoCoupling")
    plt.ylabel("Kth Neighbor EvoCoupling")
    plt.title(f"{tumor}, corr: {corr_kth}")
    plt.savefig("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/plot/Figure5_S5/plot/kth"+str(tumor)+".png")
    plt.savefig("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/plot/Figure5_S5/plot/kth"+str(tumor)+".pdf")
    
    if tumor not in corr_dict:
        corr_dict[tumor] = {}
        corr_dict[tumor]['corr_edit'] = corr_edit
        corr_dict[tumor]['corr_kth'] = corr_kth


df_corr = pd.DataFrame(columns=['Tumor', 'corr_edit', 'corr_kth'])
for tumor, corr in corr_dict.items():
    row = {'Tumor': tumor, 'corr_edit': corr['corr_edit'], 'corr_kth': corr['corr_kth']}
    df_corr = df_corr.append(row, ignore_index=True)

df_corr.to_csv('/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/plot/Figure5_S5/plot/corr_dict.csv', index=False)