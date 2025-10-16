


rm(list=ls())

library(this.path)
library(DESeq2)
library(readxl)
library(data.table)
library(dplyr)
library('ggplot2')
library(stringr)
library(tidyverse)
library(pheatmap)  # 用于作热图的包
library(ggplot2) 
library(this.path)
cur_dir2 = dirname(this.path())
setwd(cur_dir2)
getwd()
library(GEOquery)
# 加载包
library(DESeq2)
library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#GSE124305（postNAC)，GSE48075和GSE32894（Pre-NAC)

# 下载 GSE124305 数据
gse <- "GSE124305"
if(!file.exists(paste0(gse,"_eSet.Rdata"))) {
  GEO_file <- getGEO(gse,  # 需要下载的series
                     destdir = ".",  # 文件下载位置，"."表示当前目录
                     getGPL = T)  # 是否下载GPL注释文件
  save(GEO_file, file = paste0(gse,"_eSet.Rdata"))  # 将下载下来的文件保存为R语言可以处理的格式
}

gse <- "GSE48075"
if(!file.exists(paste0(gse,"_eSet.Rdata"))) {
  GEO_file <- getGEO(gse,  # 需要下载的series
                     destdir = ".",  # 文件下载位置，"."表示当前目录
                     getGPL = T)  # 是否下载GPL注释文件
  save(GEO_file, file = paste0(gse,"_eSet.Rdata"))  # 将下载下来的文件保存为R语言可以处理的格式
}
gse <- "GSE32894"
if(!file.exists(paste0(gse,"_eSet.Rdata"))) {
  GEO_file <- getGEO(gse,  # 需要下载的series
                     destdir = ".",  # 文件下载位置，"."表示当前目录
                     getGPL = T)  # 是否下载GPL注释文件
  save(GEO_file, file = paste0(gse,"_eSet.Rdata"))  # 将下载下来的文件保存为R语言可以处理的格式
}
#############################################################################################################################
加载
load("GSE124305_eSet.Rdata")
GEO_file[[1]]  # 提取GEO_file中第一个数据，有的数据有两个平台测量的数据，会有[[]]
gse124305_exp <- exprs(GEO_file[[1]])  # 提取数据中的样本基因表达矩阵
gse124305_plate <- fData(GEO_file[[1]])  #提取数据中的平台信息
gse124305_clinical <- pData(GEO_file[[1]])  # 提取数据中的样本临床信息（比如：年龄、性别、是否存活）


load("GSE48075_eSet.Rdata")
GEO_file[[1]]  # 提取GEO_file中第一个数据，有的数据有两个平台测量的数据，会有[[]]
gse48075_exp <- exprs(GEO_file[[1]])  # 提取数据中的样本基因表达矩阵
gse48075_plate <- fData(GEO_file[[1]])  #提取数据中的平台信息
gse48075_clinical <- pData(GEO_file[[1]])  # 提取数据中的样本临床信息（比如：年龄、性别、是否存活）

load("GSE32894_eSet.Rdata")
GEO_file[[1]]  # 提取GEO_file中第一个数据，有的数据有两个平台测量的数据，会有[[]]
gse32894_exp <- exprs(GEO_file[[1]])  # 提取数据中的样本基因表达矩阵
gse32894_plate <- fData(GEO_file[[1]])  #提取数据中的平台信息
gse32894_clinical <- pData(GEO_file[[1]])  # 提取数据中的样本临床信息（比如：年龄、性别、是否存活）

##############################################################
#利用下载到的平台数据进行注释

ID <- data.frame(ID_REF = gse48075_plate$ID, Gene_Symbol = gse48075_plate$ILMN_Gene)
any(duplicated(ID$ID_REF))
any(duplicated(ID$Gene_Symbol))
#ID：出现了“一个探针对应多个基因”、“探针不对应基因”、“多个探针对应一个基因”的情况
# 提取多个基因中的第一个
exp <- as.data.frame(gse48075_exp)
exp$ID_REF <- rownames(exp)  # 将探针ID添加到表达矩阵的新的一列
exp <- merge(exp, ID, by = "ID_REF") # merge函数将exp的探针id与芯片平台的探针id相匹配。merge()将二者的共同列置于第一列，其余列依次置于之后。
exp[, grep("Gene_Symbol", colnames(exp))] <- trimws(exp[, grep("Gene_Symbol", colnames(exp))])  # 去除数据头尾空格。grep("Gene_Symbol", colnames(exp))用于查找Gene_Symbol在exp的第几列。trimws用于从字符串中删除首尾空格。
#exp$Gene_Symbol <- trimws(exp$Gene_Symbol)  # 功能同上，去除数据头尾空格

exp <- na.omit(exp)  # 删除Gene_Symbol缺失的数据，解决探针不对应基因的问题
exp <- as.data.frame(exp)
write.csv(exp, "exp.gse48075.csv")
exp <- read.csv("exp.gse48075.csv", row.names = 1)

table(duplicated(exp$Gene_Symbol))  # 看一下有多少重复，TRUE为重复数。duplicated()第二次出现的基因返回的结果为FALSE。
#根据数据框某一列去重，其他的列取平均值
exp <- exp[ , -which(names(exp) %in% c("ID_REF"))]
df_new2 <- exp %>%
  group_by(Gene_Symbol) %>%  # 按基因名分组
  summarize(across(everything(), \(x) mean(x, na.rm = TRUE)))
df_new2 <- as.data.frame(df_new2)
rownames(df_new2) <- df_new2$Gene_Symbol
df_new2 <- df_new2[ , -which(names(df_new2) %in% c("Gene_Symbol"))]
exp_gse48075 <- df_new2

#################################################################################


exp <- na.omit(gse124305_exp)  # 删除Gene_Symbol缺失的数据，解决探针不对应基因的问题
exp <- as.data.frame(exp)
write.csv(exp, "exp.gse124305.csv")

exp <- read.csv("exp.gse124305.csv", row.names = 1)
   
exp_gse124305 <- exp

#######################################################################################################
ID <- data.frame(ID_REF = gse32894_plate$ID, Gene_Symbol = gse32894_plate$ILMN_Gene)
any(duplicated(ID$ID_REF))
any(duplicated(ID$Gene_Symbol))
#ID：出现了“一个探针对应多个基因”、“探针不对应基因”、“多个探针对应一个基因”的情况
# 提取多个基因中的第一个
exp <- as.data.frame(gse32894_exp)
exp$ID_REF <- rownames(exp)  # 将探针ID添加到表达矩阵的新的一列
exp <- merge(exp, ID, by = "ID_REF") # merge函数将exp的探针id与芯片平台的探针id相匹配。merge()将二者的共同列置于第一列，其余列依次置于之后。
exp[, grep("Gene_Symbol", colnames(exp))] <- trimws(exp[, grep("Gene_Symbol", colnames(exp))])  # 去除数据头尾空格。grep("Gene_Symbol", colnames(exp))用于查找Gene_Symbol在exp的第几列。trimws用于从字符串中删除首尾空格。
#exp$Gene_Symbol <- trimws(exp$Gene_Symbol)  # 功能同上，去除数据头尾空格

exp <- na.omit(exp)  # 删除Gene_Symbol缺失的数据，解决探针不对应基因的问题
exp <- as.data.frame(exp)
write.csv(exp, "exp.gse32894.csv")

exp <- read.csv("exp.gse32894.csv", row.names = 1)
table(duplicated(exp$Gene_Symbol))  # 看一下有多少重复，TRUE为重复数。duplicated()第二次出现的基因返回的结果为FALSE。
#根据数据框某一列去重，其他的列取平均值
exp <- exp[ , -which(names(exp) %in% c("ID_REF"))]
df_new2 <- exp %>%
  group_by(Gene_Symbol) %>%  # 按基因名分组
  summarize(across(everything(), \(x) mean(x, na.rm = TRUE)))
df_new2 <- as.data.frame(df_new2)
rownames(df_new2) <- df_new2$Gene_Symbol
df_new2 <- df_new2[ , -which(names(df_new2) %in% c("Gene_Symbol"))]

exp_gse32894 <- df_new2


#####################################################################################
#重复跑
library(edgeR)
library(limma)

packageVersion("edgeR")
packageVersion("limma")
#########################################################################################
exp <- read.csv("exp.gse48075.csv", row.names = 1)
table(duplicated(exp$Gene_Symbol))  # 看一下有多少重复，TRUE为重复数。duplicated()第二次出现的基因返回的结果为FALSE。
#根据数据框某一列去重，其他的列取平均值
exp <- exp[ , -which(names(exp) %in% c("ID_REF"))]
df_new2 <- exp %>%
  group_by(Gene_Symbol) %>%  # 按基因名分组
  dplyr::summarize(across(everything(), \(x) mean(x, na.rm = TRUE)))
df_new2 <- as.data.frame(df_new2)
rownames(df_new2) <- df_new2$Gene_Symbol
df_new2 <- df_new2[ , -which(names(df_new2) %in% c("Gene_Symbol"))]
exp_gse48075 <- df_new2

exp <- read.csv("exp.gse124305.csv", row.names = 1)
exp_gse124305 <- exp

exp <- read.csv("exp.gse32894.csv", row.names = 1)
table(duplicated(exp$Gene_Symbol))  # 看一下有多少重复，TRUE为重复数。duplicated()第二次出现的基因返回的结果为FALSE。
#根据数据框某一  列去重，其他的列取平均值
exp <- exp[ , -which(names(exp) %in% c("ID_REF"))]
df_new2 <- exp %>%
  group_by(Gene_Symbol) %>%  # 按基因名分组
  dplyr::summarize(across(everything(), \(x) mean(x, na.rm = TRUE)))
df_new2 <- as.data.frame(df_new2)
rownames(df_new2) <- df_new2$Gene_Symbol
df_new2 <- df_new2[ , -which(names(df_new2) %in% c("Gene_Symbol"))]
exp_gse32894 <- df_new2
# 假设 gse48075_matrix 是你的 GSE48075 表达矩阵
exp_gse48075_scaled <- scale(exp_gse48075)


common_genes <- intersect(rownames(exp_gse48075_scaled), rownames(exp_gse32894))
common_genes <- intersect(common_genes, rownames(exp_gse124305))

exp_gse124305_scaled <- exp_gse124305[common_genes, ]
exp_gse48075_scaled <- exp_gse48075_scaled[common_genes, ]
exp_gse32894_scaled <- exp_gse32894[common_genes, ] 


combined_matrix <- cbind(exp_gse124305_scaled, exp_gse48075_scaled, exp_gse32894_scaled)
write.csv(combined_matrix, "combined_matrix.csv")


batch <- c(rep("GSE124305", ncol(exp_gse124305_scaled)), 
           rep("GSE48075", ncol(exp_gse48075_scaled)), 
           rep("GSE32894", ncol(exp_gse32894_scaled)))

# 创建组别信息
group <- c(rep("postNAC", ncol(exp_gse124305_scaled)), 
           rep("preNAC", ncol(exp_gse48075_scaled) + ncol(exp_gse32894_scaled)))

#去除批次效应


combined_matrix_corrected <- removeBatchEffect(combined_matrix, batch = batch, design = model.matrix(~group))

# 创建设计矩阵
design <- model.matrix(~0 + factor(group))
colnames(design) <- c("postNAC", "preNAC")

# 适合线性模型
fit <- lmFit(combined_matrix_corrected, design)



# 创建对比矩阵
contrast_matrix <- makeContrasts(postNAC - preNAC, levels = design)

# 差异分析
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 获取差异表达基因
deg_results <- topTable(fit2, adjust.method = "BH", number = Inf)
write.csv(deg_results,file = "deg_results.csv",row.names = T)
# 筛选显著差异基因
significant_genes <- deg_results[deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 2, ]

write.csv(significant_genes,file = "significant_genes.csv",row.names = T)
######################################################################################################################################################################
##########################################################################################
#富集分析
library(ggsci) # 配色
library(ggplot2)
library(egg)

library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)
library(org.Mm.eg.db)

diff_gene <- significant_genes
gene.diff_name<-bitr(rownames(diff_gene), fromType = "SYMBOL", 
                     toType = c("SYMBOL","ENTREZID"),
                     OrgDb = "org.Hs.eg.db")

gene <- gene.diff_name$ENTREZID
ego_ALL <- enrichGO(gene = gene,#我们上面定义了
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.01,#P值可以取0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)

ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

#4、将结果保存到当前路径
ego_ALL <- as.data.frame(ego_ALL)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)
ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)#或者这样也能得到ego_ALL一样的结果
write.csv(ego_ALL,file = "ego_ALL.csv",row.names = T)
write.csv(ego_result_BP,file = "ego_result_BP.csv",row.names = T)
write.csv(ego_result_CC,file = "ego_result_CC.csv",row.names = T)
write.csv(ego_result_MF,file = "ego_result_MF.csv",row.names = T)
write.csv(ego,file = "ego.csv",row.names = T)
#################################################
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)
ego_result_CC <- ego_result_CC[order(ego_result_CC$Count,decreasing = TRUE), ]
ego_result_MF <- ego_result_MF[order(ego_result_MF$Count,decreasing = TRUE), ]
ego_result_BP <- ego_result_BP[order(ego_result_BP$Count,decreasing = TRUE), ]

#5、但是有时候我们富集到的通路太多了，不方便绘制图片，可以选择部分绘制，这时候跳过第4步，来到第5步
display_number = c(20, 20, 20)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- ego_result_BP[1:display_number[1], ]
ego_result_CC <- ego_result_CC[1:display_number[2], ]
ego_result_MF <- ego_result_MF[1:display_number[3], ]

ego_result_CC <- ego_result_CC[order(ego_result_CC$Count), ]
ego_result_MF <- ego_result_MF[order(ego_result_MF$Count), ]
ego_result_BP <- ego_result_BP[order(ego_result_BP$Count), ]

##将以上我们摘取的部分通路重新组合成数据框
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                     
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}
go_enrich_df$Description
##开始绘制GO柱状图
###横着的柱状图
go_enrich_df$Description <- factor(go_enrich_df$Description, levels = go_enrich_df$Description)
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色

p <- ggplot(data=go_enrich_df, aes(x=Description,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()+
  theme_test()+ 
  scale_y_continuous(limits = c(0,max(go_enrich_df$GeneNumber)+1),expand = c(0,0)) +
  theme(
    legend.title = element_text(color="black", size=10,),
    legend.text = element_text(color="black", size = 10),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    plot.title = element_text(size=20,colour='black'),
    axis.title = element_text(size=20,colour='black'),
    axis.text.x = element_text(angle = 0,size = 7,colour='black'),
    axis.text.y = element_text(angle = 0,size = 10,colour = 'black')
  )
p
p_modified <- egg::set_panel_size(p, width = unit(8, "in"), height = unit(7, "in"))

ggsave(
  filename = "go_all.pdf",
  plot = p_modified,
  width = 13,
  height = 10,
  units = 'in',
  dpi = 300
)
##################################


kk <- enrichKEGG(gene = gene,keyType = "kegg",organism= "hsa", qvalueCutoff = 0.05, pvalueCutoff=0.01)

#2、可视化
###柱状图
hh <- as.data.frame(kk)#
write.csv(hh,file = "kegg_result.csv",row.names = T)
rownames(hh) <- 1:nrow(hh)

hh <- hh[order(hh$Count,decreasing = TRUE), ]
hh <- hh[1:22, ]
hh <- hh[order(hh$Count), ]

##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
hh <- separate(hh, Description, into = c("Description", "Species"), sep = " - ")

hh$Description <- factor(hh$Description, levels = hh$Description)


p <- ggplot(hh,aes(x=reorder(Description, Count),y=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "#fee6ce",high ="#e6550d" )+#颜色自己可以换
  labs(title = "KEGG Pathways Enrichment",
       y = "Gene numbers", 
       x = "")+
  theme_test()+ 
  scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
  scale_y_continuous(limits = c(0,max(hh$Count)+1),expand = c(0,0)) +
  theme(
    legend.title = element_text(color="black", size=15,),
    legend.text = element_text(color="black", size = 15),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    plot.title = element_text(size=20,colour='black'),
    axis.title = element_text(size=20,colour='black'),
    axis.text.x = element_text(angle = 0,size = 8,colour='black'),
    axis.text.y = element_text(angle = 0,size = 10,colour = 'black')
  )
p
p_modified <- egg::set_panel_size(p, width = unit(8, "in"), height = unit(7, "in"))

ggsave(
  filename = "kegg.pdf",
  plot = p_modified,
  width = 12,
  height = 10,
  units = 'in',
  dpi = 300
)

##############################################################################
