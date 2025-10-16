
### cell type annotation

celltypeAnno <- function(obj, obj.markers, mouse=FALSE)
# obj: seurat3 obj, obj.markers: seurat3 "FindAllMarkers" output
# return a named vector of celltype annotations
{
  if(mouse==FALSE)
  {
    # read cell type marker list
    cellTypeMarks_csv= read.csv( file = "collectedBiomarkers_3_celltype.csv",stringsAsFactors=FALSE)
    temp_geneSymbols<- strsplit(cellTypeMarks_csv$geneSymbol, ", ")
    # create a list object to store the markers
    markers=list()
    for(i in 1:nrow(cellTypeMarks_csv))
    {
      temp_celltype=cellTypeMarks_csv$cellName[i]
      temp_markers=temp_geneSymbols[[i]]
      if(temp_celltype %in% names(markers))
      {
        markers[[temp_celltype]]=union(temp_markers,markers[[temp_celltype]])
      }else 
      {
        markers[[temp_celltype]]=c(temp_markers)
      }
    }
  }else
  {
    # celltype.markerlist2mouse() # refresh mouse celltype markers
    markers=readRDS("mouse.celltype.markers.v2.rds")
  }




  # intersect the celltype markers and the DEGs from gene expresion data
  obj.markers=obj.markers[abs(obj.markers$avg_log2FC)>log2(1.5) & obj.markers$p_val_adj<0.05, ]
  unique_markers=unique(obj.markers$gene)



  for(i in 1:(length(markers)))
  {
    markers[[i]]=intersect(markers[[i]],unique_markers)
  }
  for(i in (length(markers)):length(markers))# only keep the markers existing in expression matrix
  {
    markers[[i]]=intersect(markers[[i]],rownames(obj@assays$RNA@scale.data))
  }
  markers_union=c()
  for(i in 1:(length(markers)))
  {
    markers_union=union(markers[[i]],markers_union)
  }
  # create z-scored expression matrix for all markers
  markers_matrix_bak<-obj@assays$RNA@data[markers_union[markers_union%in%rownames(obj@assays$RNA@data)], ]
  markers_matrix<-t(scale(t(markers_matrix_bak)))
  # calculate a score for each cluster with each cell type marker
  table(colnames(markers_matrix)==names(obj@active.ident))

  clusterIDS=levels(obj@meta.data$fig.cluster)
  # clusterIDS=clusterIDS[1:length(clusterIDS)]
  # clusterScore(clusterID, celltype, obj, markers_matrix, markers)
  ### part of function "celltypeAnno"
  clusterScore<- function(clusterID, celltype ,obj, markers_matrix, markers)
  {
    markers_matrix[,obj@active.ident==clusterID]->temp_matrix
    temp_markers=markers[[celltype]]
    temp_matrix=temp_matrix[rownames(temp_matrix) %in% temp_markers,]
    if(class(temp_matrix)[1]=="matrix")
    {
      score=sum(temp_matrix)/dim(temp_matrix)[1]/dim(temp_matrix)[2]
    }else if(class(temp_matrix)[1]=="numeric")
    {
      score=sum(temp_matrix)/length(temp_matrix)
    }
    return(score)
  }


  celltypeScore<-matrix(0L, nrow = length(markers), ncol= length(clusterIDS))
  rownames(celltypeScore)=names(markers)
  colnames(celltypeScore)=as.character(clusterIDS)
  for(i in 1:length(markers))
  {
    for(j in 1:length(clusterIDS))
    {
      celltype=names(markers)[i]
      clusterID=clusterIDS[j]
      temp_score=clusterScore(clusterID, celltype, obj, markers_matrix, markers)
      celltypeScore[i,j]=temp_score
    }
  }
  highestScore=apply(celltypeScore,2,function(x) max(x, na.rm = TRUE))[]
  secondScore=apply(celltypeScore,2,function(x) max( x[x!=max(x,  na.rm = TRUE)], na.rm = TRUE ))
  judge2score=cbind(highestScore,secondScore)
  # how much the max is bigger than second
  judge2score2=apply(judge2score,1,function(x) (x[1]-x[2])/abs(x[1]))
  clusterIDS_celltype=clusterIDS
  # if second cell type's score is >20% lower than best's, only show best's score, else show both scores
  for(i in 1:length(clusterIDS_celltype))
  {
    temp_col=celltypeScore[,i]
    bestCelltype_score=temp_col[temp_col==highestScore[i]]
    bestCelltype_score=bestCelltype_score[!is.na(bestCelltype_score)]
    bestCelltype=names(bestCelltype_score)
    secondCelltype_score=temp_col[temp_col==secondScore[i]]
    secondCelltype_score=secondCelltype_score[!is.na(secondCelltype_score)]
    secondCelltype=names(secondCelltype_score)

    temp_name=paste(bestCelltype,": ",round(bestCelltype_score,2),sep="")
    temp_name=paste(temp_name, collapse=" ")
    if(judge2score2[i]<0.1)
    {
      temp_name_second=paste(secondCelltype,": ",round(secondCelltype_score,2),sep="")
      temp_name_second=paste(temp_name_second, collapse=" ")
      temp_name=paste(temp_name, " ", temp_name_second, sep="")
    }
    clusterIDS_celltype[i]=temp_name
  }
  print("Cell type score matrix:") 
  print(celltypeScore)
  print("Decided cell type:")
  print(clusterIDS_celltype)

  # check if there are multiple same name
  table(clusterIDS_celltype)
  # simplify the names of decided cell types
  rm.position=regexpr("[:_]", clusterIDS_celltype, perl=TRUE)
  clusterIDS_celltype=substring(clusterIDS_celltype, 1, (rm.position-1))
  table(clusterIDS_celltype)

  ## show the cluster's celltype in T-sne map
  current.cluster.ids <- clusterIDS
  new.cluster.ids <- clusterIDS_celltype

  obj@meta.data$celltype <- plyr::mapvalues(x = obj@meta.data$fig.cluster, from = current.cluster.ids, to = new.cluster.ids)
  temp.type=as.character(obj@meta.data$celltype)
  names(temp.type)=rownames(obj@meta.data)
  return(temp.type)
}

