require(MetaVolcanoR)
.libPaths(c("/home/DewangZhou/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))
library(clusterProfiler)
#BiocManager::install("MetaVolcanoR",lib="/home/yanzeqin/R/x86_64-pc-linux-gnu-library/4.0")
#BiocManager::install("clusterProfiler",lib="/home/yanzeqin/R/x86_64-pc-linux-gnu-library/4.0")
require(ggplot2)
require(ggrepel)
require(biomaRt)
require(msigdbr)
library(enrichplot)
require(cowplot)
require(dplyr)
require(SuperExactTest)

NUM_CELLS_THRESH = 10
PERCENT_UNIQUE_THRESH = 0.05
PERCENT_UNSATURATED_TARGETS_THRESH = 0.2

data_directory = "/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/"
tumor_clone_statistics = read.table(paste0(data_directory, "/allele/tumor_statistics.tsv"), sep='\t', header=T, row.names=1)
#这一步省略
filt.tumors = rownames(subset(tumor_clone_statistics, NumCells < NUM_CELLS_THRESH | 
                                PercentUnique <= PERCENT_UNIQUE_THRESH |
                                PercentUnsaturatedTargets <= PERCENT_UNSATURATED_TARGETS_THRESH))

tumor_list <- read.table(paste0(data_directory,"/allele/tumor_list.txt"), sep='\t', header=F, row.names=1, stringsAsFactors = F)$Tumors
a <- read.table(paste0(data_directory,"/allele/tumor_list.txt"), sep='\t',header= T, stringsAsFactors = F)
tumor_list = a$Tumor



de.results = list()

beta.df = data.frame()

home.dir = paste0(data_directory, "/allele/trees")
for (tumor in tumor_list) {
  
  if (grepl('MGH', tumor, fixed=T)) {
    
    gene.fp = paste(home.dir, paste0("linregress.", tumor, ".log2fc.txt"), sep="/")
    
    if (file.exists(gene.fp)) {
      
      res = read.table(gene.fp, sep='\t', header=T)
      
      if (!'genes' %in% colnames(res)) {
        res$genes = res$X
      }
      
      de.results[[tumor]] = res
    }
  }
}
##################################################################################################
#火山图
tumor = 'CTR_MGH_0'
fp = paste0(home.dir, "/linregress.", tumor, ".log2fc.txt")
LOG2FC = log2(1.5)
FDR_THRESH = 0.05

savefig = FALSE

DISC_THRESH = 50

genes.df = read.table(fp, sep='\t', header=T)
genes.df$log10qval = -log10(genes.df$FDR)

# remove genes whose pvalue is 0
genes.df[is.infinite(genes.df$log10qval),"log10qval"] = 300
genes.df$log10qval = unlist(lapply(genes.df$log10qval, function(x) min(300, x)))

genes.df$disc = apply(genes.df, 1, function(x) as.numeric(x["log10qval"]) * abs(as.numeric(x["log2fc"])))                              


hits = genes.df[genes.df$disc >= DISC_THRESH, ]
hit.genes = hits$genes

pos = hits[hits$log2fc >= 0, ]
neg = hits[hits$log2fc < 0, ]

pos.ordered = pos[order(-pos$disc), ]
neg.ordered = neg[order(-neg$disc), ]


pos.genes = as.character(pos.ordered$genes)[1:20]
neg.genes = as.character(neg.ordered$genes)[1:20]

print(dim(subset(genes.df, log2fc > LOG2FC & FDR < FDR_THRESH)))
print(dim(subset(genes.df, log2fc < -1*LOG2FC & FDR < FDR_THRESH)))


ggplot(genes.df, aes(log2fc, log10qval, label=genes)) + 
  geom_point(data = subset(genes.df, abs(log2fc) < LOG2FC | FDR >= FDR_THRESH), color='black') +
  geom_point(data = subset(genes.df, log2fc > LOG2FC & FDR < FDR_THRESH), color='red') + 
  geom_point(data = subset(genes.df, log2fc < -1*LOG2FC & FDR < FDR_THRESH), color='blue') +
  geom_text_repel(data = subset(genes.df, genes %in% hit.genes), size = 2.5) + 
  labs(x = 'Log2FC', y = 'Log10(FDR)', title = paste0('Regression Volcano, ', tumor)) +
  theme_classic() + theme(aspect.ratio = 1)

########################################################################################################
#Majority Vote Analysis
meta_degs_vc <- votecount_mv(diffexp=de.results,
                             pcriteria='FDR', 
                             foldchangecol='log2fc',
                             genenamecol='genes',
                             pvalue=0.05,
                             foldchange=log2(1.5),
                             metathr=0.01,
                             collaps=T,
                             jobname='MetaVolcano',
                             outputfolder='.',
                             draw='PDF')

res = meta_degs_vc@metaresult

head(res)
SIGN_CONSISTENCY = 2
SIGN_TO_FREQ = 0.5
FREQ = 0

res = meta_degs_vc@metaresult

res$ddeg = as.numeric(res$ddeg)
res$ndeg = as.numeric(res$ndeg)
res$sign_to_freq = apply(res, 1, function(x) abs(as.numeric(x[['ddeg']])) / max(1, as.numeric(x['ndeg'])) )

sig.genes = res[which(res$sign_to_freq > SIGN_TO_FREQ),]
up.genes = as.character(sig.genes[which(sig.genes$ddeg >= SIGN_CONSISTENCY),"genes"])
down.genes = as.character(sig.genes[which(sig.genes$ddeg <= (-1*SIGN_CONSISTENCY)),"genes"])

hit.genes = c(up.genes, down.genes)

rownames(res) = as.character(res$genes)

res$degvcount <- as.character(res$degvcount)
res$degvcount <- '1.Unperturbed'
res[up.genes, 'degvcount'] <-'2.Up'
res[down.genes, 'degvcount'] <- '0.Down'
res$degvcount = factor(res$degvcount, levels=c("0.Down", "1.Unperturbed", "2.Up"))


g = ggplot(res, aes(x = ddeg, y = ndeg)) +
  geom_jitter(aes(color = degvcount), 
              cex=0.5, width=0.45, height=0.45) +
  labs(x = "Sign Consistency",
       y = "Number of Studies as Differentially Expressed",
       title = "Majority Vote Differential Expression Consensus, sgNT")

g + theme_classic() +
  theme(panel.border= element_blank()) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) +
  theme(axis.line.x = element_line(color = "black", size = 0.6, 
                                   lineend = "square"),
        axis.line.y = element_line(color = "black", size = 0.6, 
                                   lineend = "square")) +
  theme(legend.position = "none") +
  scale_color_manual(values=c("#377EB8", "grey", "#E41A1C"))

up.genes.df = res[up.genes, c("genes", 'ndeg', 'ddeg', 'sign_to_freq')]
colnames(up.genes.df) <- c("genes", "FreqDE", "ConsistencyDE", "ConsistencyToFreq")

down.genes.df = res[down.genes, c("genes", 'ndeg', 'ddeg', 'sign_to_freq')]
colnames(down.genes.df) <- c("genes", "FreqDE", "ConsistencyDE", "ConsistencyToFreq")


print(paste0('Found ', nrow(up.genes.df), ' up-regulated genes and ', nrow(down.genes.df), ' down-regulated genes.'))


write.table(down.genes.df,"/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/plot/Figure3_S3/data/majority_vote.down_genes.sgNT.txt",  sep="\t", row.names=F,quote=FALSE)
write.table(up.genes.df,"/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/plot/Figure3_S3/data/majority_vote.up_genes.sgNT.txt",  sep="\t", row.names=F,quote=FALSE)

##############################################################################################################
#GO Analysis of up-regulated genes

m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
m_bp = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
bp.up = enricher(gene = up.genes, TERM2GENE = m_bp)
#挑选基因通路
#up.programs = c("GO_RIBOSOME_BIOGENESIS",
#                "GO_TRANSLATIONAL_INITIATION",
#                "GO_MRNA_PROCESSING",
#                "GO_STEM_CELL_DIFFERENTIATION", 
#                "GO_REGULATION_OF_RESPONSE_TO_WOUNDING", 
#                "GO_AMEBOIDAL_TYPE_CELL_MIGRATION",
#                "GO_RIBOSE_PHOSPHATE_METABOLIC_PROCESS",
#                "GO_RNA_DEPENDENT_DNA_BIOSYNTHETIC_PROCESS",
#                "GO_OSSIFICATION",
#                "GO_GLIAL_CELL_DEVELOPMENT",
#                "GO_EPITHELIAL_CELL_PROLIFERATION",
#                "GO_LAMELLIPODIUM_ORGANIZATION",
#                "GO_GLIAL_CELL_DIFFERENTIATION",
#                "GO_SKIN_DEVELOPMENT")
#results = bp.up@result[up.programs, ]
results = bp.up@result
results
idx <- order(results[['p.adjust']], decreasing = F)
results$Description <- factor(results$Description,
                              levels=rev(unique(results$Description[idx])))

ggplot(results, aes_string(x='p.adjust', y="Description", size='Count')) +
  geom_point() +
  guides(size  = guide_legend(order = 1)) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 7))

##############################################################################################################################################################################################
meta = read.table("./data/marjanovic_meta.tsv", sep='\t', row.names=1, header=T)
marjanovic_sigscores = read.table('./data/marjanovic_signatures.tsv', sep='\t', row.names=1, header=T)

meta$FitnessSignature = marjanovic_sigscores[rownames(meta), 'FitnessSignature']

mean_fitnesses = meta %>%
  group_by(timesimple) %>%
  summarize(mean_fitness = mean(FitnessSignature))

mean_fitness_kp = mean_fitnesses[mean_fitnesses$timesimple %in% c('06_KP_12w_ND', '07_KP_20w_ND', '08_KP_30w_ND'),]


ggplot(mean_fitness_kp, aes(x=timesimple, y = mean_fitness))+
  geom_col(fill='black') + 
  theme_bw() +
  scale_x_discrete(breaks=c("06_KP_12w_ND","07_KP_20w_ND","08_KP_30w_ND"),
                   labels=c('12 weeks', '20 weeks', '30 weeks')) + 
  labs(y = 'Mean FitnessSignature', x = 'Timepoint', title = 'Mean FitnessSignature of KP Timepoints from Marjanovic et al (2020)')

sessionInfo()
