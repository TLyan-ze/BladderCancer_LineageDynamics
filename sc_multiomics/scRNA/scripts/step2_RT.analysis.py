import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["font.sans-serif"] = "Arial"

sc.settings.set_figure_params(dpi=50, dpi_save=300, figsize=(4, 4))


data_directory = "/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/"

adata = sc.read_h5ad(f'/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/expression/RT.h5ad')

adata.obs['A'] =adata.obs.index
adata.obs['B'] = adata.obs['A'].str.split('-',expand=True)[0]
adata.obs.set_index('B',inplace=True)
adata.obs.index.name = None
adata.obs['A'] =adata.obs.index
print(adata.obs)

allele_table_unfiltered = pd.read_csv(f"/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele_RT/01_metainfo/AP.alleleTable.unfiltered.txt2", sep='\t')
allele_table_unfiltered

new_df = allele_table_unfiltered.loc[:, ["cellBC","Tumor"]]
df = new_df.drop_duplicates(subset=['cellBC'])
a=adata.obs
merged_df = pd.merge(a,df, left_on='A', right_on='cellBC', how='left')
merged_df = merged_df.fillna('NA')
merged_df.set_index('A',inplace=True)
merged_df.index.name = None
merged_df['A'] =merged_df.index
adata.obs=merged_df
print(adata.obs)
##################

#### 简单的过滤一下细胞和基因，一个细胞至少表达300个基因，一个基因至少在10个细胞中表达。
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(adata)
#计算线粒体基因的比例，线粒体基因以MT开头
adata.var['mt'] = adata.var_names.str.startswith(r'MT-',r'mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
#可视化count数量、基因数量、线粒体基因百分比
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
plt.savefig('/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/expression/plot_RT/01-qc_violin1.png')
sc.settings.set_figure_params(dpi_save=300)
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8, 3))
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=ax[0], show=False)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=ax[1], show=False)
plt.subplots_adjust(wspace=.4)


upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
upper_lim


adata = adata[adata.obs.n_genes_by_counts < 11000]
adata = adata[adata.obs.n_genes_by_counts > 3000]


adata = adata[adata.obs.pct_counts_mt < 20]

#adata = adata[adata.obs.pct_counts_ribo < 35]
adata

#####################################################################
print(1)
#可能有多种原因使得每个细胞在测序时候总的count不同，为了细胞间可以比较因此将他们都统一为10000。
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#一般筛选1000-5000个高变基因（HVG），可以降低数据维度。
#HVG用于降维、聚类及其可视化等。
#同时也在adata存一份全部基因的rawdata用于marker鉴定、差异分析、细胞周期分数计算、基因表达量可视化等。

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
sc.pl.highly_variable_genes(adata)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]

#regress_out用普通的线性回归消除pct_counts_Mito变量的影响。
#scale可选，把数据中心化
sc.pp.regress_out(adata, ['pct_counts_mt',])
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, use_highly_variable=True, svd_solver='arpack',n_comps=50)
sc.pl.pca(adata, color="sample",legend_loc="on data",legend_fontsize="small")
sc.pl.pca_variance_ratio(adata,n_pcs=50)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.pl.umap(adata, color="sample",legend_loc="on data",legend_fontsize="small")

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata)
sc.pl.umap(adata, color='leiden',)

########################
print(2)
import scanpy as sc
#使用leiden_sub聚类
sc.tl.leiden(adata, resolution=0.7, key_added='leiden_sub', random_state=0)
#使用 UMAP 算法将聚类结果可视化，并使用颜色编码表示不同的聚类
sc.pl.umap(adata, color='leiden_sub')

sc.pl.umap(adata, color='sample')
sc.tl.tsne(adata)
sc.pl.tsne(adata, color='sample')

sc.pl.umap(adata, color="sample",legend_loc="on data",legend_fontsize="small")
plt.savefig('/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/expression/plot_RT/03-clsters_violin1.png')
sc.pl.umap(adata, color='leiden_sub',legend_loc="on data",legend_fontsize="small")
plt.savefig('/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/expression/plot_RT/03-clsters_violin2.png')
adata.var_names_make_unique()
adata.obs_names_make_unique()
adata.write('scRNA4_2000_RT.h5ad')



