rm(list=ls())
library(this.path)
cur_dir2 = dirname(this.path())
setwd(cur_dir2)
getwd()

library(ggsci) # 配色
library(ggplot2)
library(egg)

library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)#画图包

all_gene <-read.csv(file="MI_IGFKO_2.csv")



############################################################################
#上调
diff_gene <- subset(all_gene,pvalue < 0.05 & abs(log2FC) > 1)

gene.diff_name<-bitr(diff_gene$SYMBOL, fromType = "SYMBOL", 
                     toType = c("SYMBOL","ENTREZID"),
                     OrgDb = "org.Hs.eg.db")

gene <- diff_gene$ENTREZID
#3、GO富集
##CC表示细胞组分，MF表示分子功能，BP表示生物学过程，ALL表示同时富集三种过程，选自己需要的,我一般是做BP,MF,CC这3组再合并成一个数据框，方便后续摘取部分通路绘图。
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


ego_result_CC <- ego_result_CC[order(ego_result_CC$Count,decreasing = TRUE), ]
ego_result_MF <- ego_result_MF[order(ego_result_MF$Count,decreasing = TRUE), ]

#
ego_result_BP <- as.data.frame(ego_BP)
ego_result_BP <- ego_result_BP[order(ego_result_BP$Count,decreasing = TRUE), ]
ego_result_BP <- ego_result_BP[1:20, ]
rownames(ego_result_BP) <- 1:nrow(ego_result_BP)
ego_result_BP$order=factor(rev(as.integer(rownames(ego_result_BP))),labels = rev(ego_result_BP$Description))

ego_result_BP$p.adjust <- as.numeric(ego_result_BP$p.adjust)
p <- ggplot(ego_result_BP,aes(x=reorder(order, Count),y=Count,fill = -log10(p.adjust)))+
  geom_bar(stat = "identity",width=0.8)+####柱子宽度
  coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "#fee6ce",high ="#e6550d" )+#颜色自己可以换
  labs(title = "Pathway Enrichment",y = "Number of Upregulated Genes", x = "",fill="-log10(p.adjust)")+
  theme_test()+ 
  scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
  scale_y_continuous(limits = c(0,max(ego_result_BP$Count)+1),expand = c(0,0)) +
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




















###########################################################################################################

#1、KEGG富集
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism= "human", qvalueCutoff = 0.05, pvalueCutoff=0.05)

#2、可视化
###柱状图
hh <- as.data.frame(kk)#
write.csv(hh,file = "kegg_result.csv",row.names = T)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))

#x = reorder(Description2, Genes_in_Overlap_k)

p <- ggplot(hh,aes(x=reorder(order, Count),y=Count,fill=p.adjust))+
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
###气泡图
hh <- as.data.frame(kk)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
ggplot(hh,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=-1*p.adjust))+# 修改点的大小
  scale_color_gradient(low="green",high = "red")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()
