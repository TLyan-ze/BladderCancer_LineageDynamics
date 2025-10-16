
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
cur_dir2 = dirname(this.path())
setwd(cur_dir2)
getwd()


#library(GenomicFeatures)
#########################################################################################
#数据框处理，列为样本，行为基因名
exp_info <- read.csv("GEMM mus_counts.tsv",sep = "\t" ) 


exp_info$Gene_Name <- sapply(strsplit(exp_info$X, "\\|"), function(x) tail(x, n=3)[1])

exp_info <- exp_info[ , -which(names(exp_info) %in% c("X"))]


#根据数据框某一列去重，其他的列取平均值

df_new2 <- exp_info %>%
  group_by(Gene_Name) %>%  # 按基因名分组
  summarize(across(everything(), \(x) mean(x, na.rm = TRUE)))
df_new2 <- as.data.frame(df_new2)
rownames(df_new2) <- df_new2$Gene_Name
df_new2 <- df_new2[ , -which(names(df_new2) %in% c("Gene_Name"))]

#去掉表达量为0
df_new2 <- df_new2[rowMeans(df_new2)>1,]    


condition <- factor(c(rep("control",3),rep("treat",3)))

colData <- data.frame(row.names=colnames(df_new2), condition)

df_new2_int <- round(df_new2)
dds <- DESeqDataSetFromMatrix(countData = df_new2_int, colData = colData, design = ~ condition)
head(dds)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
res <- results(dds1)

summary(res) 
# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

# 依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
write.csv(res1,file = "GEMM_mus_DEseq2_all.csv",row.names = T)

res1_up<- res1[which(res1$log2FoldChange >= 1 & res1$pvalue < 0.05),]      # 表达量显著上升的基因
res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$pvalue < 0.05),]    # 表达量显著下降的基因
res1_total <- rbind(res1_up,res1_down)

write.csv(res1_total,file = "GEMM_mus_DEseq2.csv",row.names = T)


library(ggrepel)
genes<- res1
# 根据上调、下调、不变为基因添加颜色信息
genes$color <- ifelse(genes$padj<0.05 & abs(genes$log2FoldChange)>= 1,ifelse(genes$log2FoldChange > 1,'red','blue'),'gray')
color <- c(red = "red",gray = "gray",blue = "blue")

p <- ggplot(
  # 指定数据、映射、颜色
  genes, aes(log2FoldChange, -log10(padj), col = color)) +  
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  # 辅助线
  labs(x="log2 (fold change)",y="-log10 (q-value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  # 图例
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)+
          # 注释
          geom_text_repel(
            data = subset(genes, padj < 1e-100 & abs(genes$log2FoldChange) >= 10),
            aes(label = rownames(genes)),
            size = 5,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines"))
  )
p


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

diff_gene <- res1_total
gene.diff_name<-bitr(rownames(diff_gene), fromType = "SYMBOL", 
                     toType = c("SYMBOL","ENTREZID"),
                     OrgDb = "org.Mm.eg.db")

gene <- gene.diff_name$ENTREZID
ego_ALL <- enrichGO(gene = gene,#我们上面定义了
                    OrgDb=org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.01,#P值可以取0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)

ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Mm.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Mm.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Mm.eg.db,
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


kk <- enrichKEGG(gene = gene,keyType = "kegg",organism= "mmu", qvalueCutoff = 0.05, pvalueCutoff=0.01)

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
########################################################################################################
#gsea
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
library(enrichplot)


DESeq2 <- res1
head(DESeq2)
length(rownames(DESeq2))
symbol<- rownames(DESeq2)
entrez<- bitr(symbol,
              
              fromType= "SYMBOL",#现有的ID类型
              
              toType= "ENTREZID",#需转换的ID类型
              
              OrgDb= "org.Mm.eg.db")
head(entrez)
genelist<- DESeq2$log2FoldChange
names(genelist) <- rownames(DESeq2)

genelist<- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
head(genelist)
genelist<- sort(genelist, decreasing = T)

head(genelist)

GO_ges<- gseGO(geneList = genelist,
               
               OrgDb= org.Mm.eg.db,
               
               ont= "ALL", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
               pAdjustMethod = "BH",
      
               pvalueCutoff= 0.05,
               )

GO_ges<- setReadable(GO_ges, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

GO_ges_result<- GO_ges@result

write.csv(GO_ges_result, file = c('GO_gsea_result2.csv'))

p4<- gseaplot2(GO_ges,
               
               geneSetID= 30,
               
               color= "red",
               
               rel_heights= c(1.5, 0.5, 1), #子图高度
               
               subplots= 1:3, #显示哪些子图
               
               pvalue_table= T, #是否显示pvalue表
               
               title= GO_ges@result$Description[30],
               
               ES_geom= "line") #"dot"将线转换为点

p4
