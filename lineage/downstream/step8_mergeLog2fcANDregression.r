


data_directory = "/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/"

a <- read.table(paste0(data_directory,"/allele/tumor_list.txt"), sep='\t',header= T, stringsAsFactors = F)
tumor_list = a$Tumor
filt.tumor <- list("CTR_MGH_0")

for (tumor in tumor_list) {
  if (grepl('MGH', tumor, fixed=T)){
    print(1)
    log2fc=paste0("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/fitnesses/fitness_log2fc.", tumor, ".txt")
    log2fc_df = read.table(fp, sep='\t', header=T)
    linregress=paste0("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/fitnesses/linregress.", tumor, ".txt")
    linregress_df = read.table(linregress, sep='\t', header=T)
    merged_trcr <- merge(log2fc_df, linregress_df, by.x="genes", by.y="genes")
    adjusted_p_values <- p.adjust(merged_trcr$pval, method = "BH")
    merged_trcr$FDR <- adjusted_p_values
    output_df =paste0("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele/trees/linregress.", tumor, ".log2fc.txt")
    write.table(merged_trcr,output_df,  sep="\t", row.names=F,quote=FALSE)
    }
}



