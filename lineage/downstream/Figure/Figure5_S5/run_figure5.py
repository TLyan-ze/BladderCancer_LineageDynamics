#coding:utf-8
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
#from cassiopeia.Analysis import small_parsimony

sys.path.append("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/KPTracer-release-main/reproducibility/Figure5_S5/scripts/")
import tree_utilities


data_directory = "/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/"

adata = sc.read_h5ad(f'{data_directory}/expression/scRNA4.h5ad')
print(adata)

tumor= 'M_2'

graph = tree_utilities.prepare_tumor_tree(tumor, adata,tree_dir = f"/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele/trees/",
                                            FILTER_PROP=0.025,
                                            column='seurat_clusters')



data_directory = "/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/KPTracer-Data"

adata2 = sc.read_h5ad(f'{data_directory}/expression/adata_processed.nt.h5ad')
#leiden_sub
print(graph)