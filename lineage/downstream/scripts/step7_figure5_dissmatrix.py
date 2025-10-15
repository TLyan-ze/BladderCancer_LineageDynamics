import os
import sys

import ete3
import networkx as nx
import numba
import numpy as np
import pandas as pd
import scanpy as sc
from tqdm.auto import tqdm
# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))

# 将当前工作目录设置为脚本所在目录
os.chdir(script_dir)

import tree_utilities




data_directory = script_dir
adata = sc.read_h5ad(
    f'/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/expression/scRNA4_3000_25.h5ad'
)
tumor_list = pd.read_csv(f"{data_directory}/tumor_list.txt", sep='\t')['Tumor'].values

FILTER_PROP = 0.025
PREPROCESS = True

cluster_column = "leiden_sub"
print(adata.obs)
uniq_states = adata.obs[cluster_column].unique()

tumor_z_df = {}
observed_df = {}

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
for tumor in tqdm(tumors):

    print(tumor)

    graph = tree_utilities.prepare_tumor_tree(
        tumor,
        adata,
        tree_dir = f"{data_directory}/trees",
        column=cluster_column,
        preprocess=PREPROCESS,
    )

    if PREPROCESS:
        for u, v in graph.edges():
            _length = graph[u][v]["length"]
            if _length > 0:
                _length = 1
            graph[u][v]["length"] = _length

    (
        phylogenetic_distance_matrix,
        leaf_pairs,
        edit_distance_matrix,
        tree_diameter,
        n_targets,
    ) = tree_utilities.compute_pairwise_dist_nx(graph)

    np.fill_diagonal(phylogenetic_distance_matrix.values, 0)
    np.fill_diagonal(edit_distance_matrix.values, 0)

    phylogenetic_distance_matrix.to_csv(
        f"{data_directory}/Figure/Figure5_S5/data/{tumor}_phylogenetic_distance_matrix.tsv",
        sep="\t",
    )

    edit_distance_matrix.to_csv(
        f"{data_directory}/Figure/Figure5_S5/data/{tumor}_edit_distance_matrix.tsv",
        sep="\t",
    )
