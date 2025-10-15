import hotspot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, cophenet, fclusterdata
from scipy.spatial.distance import squareform
import seaborn as sns



# read in adata
data_directory = "/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data"
adata = sc.read(f'{data_directory}/expression/scRNA4_3000_25.h5ad')

mv_up_genes = pd.read_csv("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele_new/Figure/Figure3_S3/data/majority_vote.up_genes.sgNT.txt", sep='\t')['genes'].values
mv_down_genes = pd.read_csv("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele_new/Figure/Figure3_S3/data/majority_vote.down_genes.sgNT.txt", sep='\t')['genes'].values
#####################################################
#处理数据
scale_factor = np.median(np.array(adata.X.sum(axis=1)))
filt = adata.obs.apply(lambda x: 'Normal' not in x.Tumor and type(x.Tumor) == str, axis=1)
adata_filt = adata[filt,:]
mv_genes = mv_up_genes
mv_genes = [gene for gene in mv_genes if gene in adata_filt.var_names]
adata_mv = adata_filt[:,mv_genes]
counts = pd.DataFrame(adata_mv.X.astype("float64"), index = adata_mv.obs_names, columns = adata_mv.var_names).T
umi_counts = counts.sum(axis=0)
latent = pd.DataFrame(adata_mv.obsm['X_pca'], index = adata_mv.obs_names)

#创建 hotspot.Hotspot 对象：基于经过筛选和处理后的数据，使用 hotspot.Hotspot 类创建 hs 对象，并提供 counts、model、latent 和 umi_counts 参数。
hs = hotspot.Hotspot(counts, model='danb', latent=latent, umi_counts=umi_counts)
hs = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model='danb',
    latent_obsm_key="X_pca",
    umi_counts_obs_key="total_counts"
)
#创建 KNN 图：使用 create_knn_graph 方法创建一个 KNN 图，通过设置参数如 weighted_graph 和 n_neighbors 来构建图。
hs.create_knn_graph(weighted_graph=False, n_neighbors=int(np.sqrt(counts.shape[1])))
#计算自相关性：使用 compute_autocorrelations 方法计算自相关性，并将结果保存在 hs_results 变量中。
hs_results = hs.compute_autocorrelations(jobs=10)
hs_results.head(15)
###################
#对自相关性结果进行分析：从 hs_results 中选择具有统计显著性的基因，并按照 FDR 和 Z 值进行排序
hs_genes = hs_results.loc[hs_results.FDR < 0.05].sort_values('Z', ascending=False).index
#进行局部相关性计算：使用 compute_local_correlations 方法计算局部相关性，并将结果保存在 lcz 变量中。
lcz = hs.compute_local_correlations(hs_genes, jobs=20)
#定义一个子聚类函数 subcluster：该函数根据 hs 对象的局部相关性矩阵对指定的基因进行子聚类分析。
def subcluster(hs, genes, desired_num_clusters = 2):
    
    dd = hs.local_correlation_z.loc[genes, genes].copy().values
    np.fill_diagonal(dd, 0)
    condensed = squareform(dd)*-1
    offset = condensed.min() * -1
    condensed += offset
    Z = linkage(condensed, method='average')
    
    clusters = fcluster(Z, t=desired_num_clusters, criterion='maxclust')
    return clusters

hs.local_correlation_z = lcz
#创建模块：通过调用 create_modules 方法，创建模块并根据一些参考值进行过滤
modules = hs.create_modules(
    min_gene_threshold=100, core_only=False, fdr_threshold=0.05
)
#进行层次聚类分析：使用 linkage 函数计算局部相关性矩阵的层次聚类，并根据阈值 NUM_CLUSTERS 成为不同的聚类。
NUM_CLUSTERS = 3
dd = hs.local_correlation_z.copy().values
np.fill_diagonal(dd, 0)
condensed = squareform(dd)*-1
offset = condensed.min() * -1
condensed += offset
Z = linkage(condensed, method='ward')
clusters = fcluster(hs.linkage, t = NUM_CLUSTERS, criterion='maxclust')
#对聚类结果进行处理：将聚类结果保存在 clusters 变量中，并将其转换为 pd.DataFrame 对象
modules = pd.DataFrame(clusters, index=hs.local_correlation_z.index).iloc[:,0]
hs.modules = None
hs.modules = modules
hs.modules.name = 'Module'
#可视化结果：使用 plot_local_correlations 方法绘制局部相关性图像，并将图像保存为 hotput.png 文件。
hs.plot_local_correlations(vmin=-15, vmax=15)
plt.show()
plt.savefig("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele_new/Figure/Figure3_S3/hotput.png")