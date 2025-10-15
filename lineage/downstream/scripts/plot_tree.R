library(tidyverse)
library(ggtree)
library(treeio)
library(ggsci)


library(ggtreeExtra)
library(ggstar)
library(ggtree)
library(ggplot2)
library(treeio)
library(ggnewscale)

#tree <- read.newick("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele/trees/Met_MGH_2_tree.nwk")

tree <- read.tree("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele/trees/Met_MGH_2_tree.nwk")


dat1 <- read.csv("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/expression/scRNA3_obs.csv",header = T,sep = "\t")
knitr::kable(head(dat1))
dat1$leiden <- as.character(dat1$leiden)


p <- ggtree(tree, layout="fan", 
            branch.length = "none",
            geom_tiplab2(size=3),
            open.angle=0, size=0.1)
p


p2 <- p + 
  new_scale_fill()+
  geom_fruit(
    data=dat1,
    geom=geom_col,
    mapping=aes(y=cellBC, x="10", fill=leiden),  #The 'Abundance' of 'dat1' will be mapped to x
  )
p2

