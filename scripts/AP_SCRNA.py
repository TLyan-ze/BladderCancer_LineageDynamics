import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib as mpl
params = {
    'font.family': "Helvetica",
    'figure.dpi': 300
   }
mpl.rcParams.update(params)
mpl.rc('savefig', dpi=300)

data_directory = "/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/KPTracer-Data"
adata = sc.read_h5ad("/WorkDir2/yanzeqin/ScRNA-seq/test.h5ad")

print(1)
print(adata.var)
print(2)
print(adata)
print(3)
print(adata.obs.shape) # n个细胞
print(4)
print(adata.var.shape) # n个基因
print(5)
print(adata.to_df().shape)
print(6)
print(adata.obs.head())
print(7)
print(adata.var.head())
#sc.pl.umap(adata, color='leiden_sub')

#sigscores = pd.read_csv(f"../Figure6_S6/data/fitness_signature_scores.tsv", sep='\t', index_col = 0, usecols=['FitnessSignature_NT'])