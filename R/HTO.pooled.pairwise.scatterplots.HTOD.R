# Plot HTO barcodes in pairwise manner, raw counts
#rm(list = ls())
library(Matrix)
library(Biobase)
#name <- "CZI_HTO_pooled"
#setwd("/Users/lawlon/Documents/CZI/CITE-seq/PBMC_2018/Pooled_Runs/HTO_Analysis/")
# read in all hto

HTO.pooled.pairwise.scatterplots.HTOD=function(mtx_file,clss_file){

hto_sel <- readRDS(mtx_file)
#bars <- c("Bcell-TGATGGCCTATTGGG", "Mono-AGTAAGTTCAGCGTA", "Tcell-TTCCGCCTCTCTTTG", "control-GTCAACTCTTTAGCG")
#hto_sel <- hto[bars, ]

#"CZI.PMBC.HTO.pooled.HTO.matrix.Rds"
    
#pr_name=strsplit(x=file_name,split=".HTO.pooled")[[1]][1] 
    
name=strsplit(x=mtx_file,split=".pooled")[[1]][1]    
    
# pairwise scatterplots
pdf(file = paste(name, "raw_counts_pairwise_scatterplots.pdf", sep = "_"), onefile = T)
for (i in 1:(nrow(hto_sel)-1)) {
  j = i + 1
  for (j in j:nrow(hto_sel)) {
    # raw
    #jpeg(file = paste(name,i,"vs",j, "raw_counts_pairwise_scatterplots.jpg", sep = "_"))
    plot(x = as.numeric(hto_sel[i,]), y = as.numeric(hto_sel[j,]),
         main = "HTO Raw Read Counts", xlab = rownames(hto_sel)[i], ylab = rownames(hto_sel)[j], cex = 0.75)
   # dev.off()
    # log scaled
    #jpeg(file = paste(name,i,"vs",j, "raw_log_counts_pairwise_scatterplots.jpg", sep = "_"))

    plot(x = log2(as.numeric(hto_sel[i,])+1), y = log2(as.numeric(hto_sel[j,])+1),
         main = "HTO Log2 Raw Read Counts", xlab = rownames(hto_sel)[i], ylab = rownames(hto_sel)[j], cex = 0.75)
   # dev.off()
  }
    

}
dev.off()

#CZI.PBMC.HTO.sample.classifications.GMM.csv

#gr_trt_fldr=paste0("/projects/ucar-lab/danaco/GroundTruths/",pr_name,"/")
    
#clss_name=paste(name,"sample.classifications.GMM.csv",sep=".")
    
# load in hto annotations
hto_anns <- read.csv(clss_file)
#table(rownames(hto_anns) == colnames(hto_sel))

pdf(file = paste(name, "labeled_pairwise_scatterplots_ggplots.pdf", sep = "_"), onefile = T)

for (i in 1:(nrow(hto_sel)-1)) {
  j = i + 1
  for (j in j:nrow(hto_sel)) {
    #jpeg(file = paste(name, i, "vs", j, "pairwise_scatterplots_ggplots.jpeg", sep = "_"))
    # create a dataframe for ggplot
    exp_df <- data.frame(Cell_Barcode = colnames(hto_sel),
                         Annotation = hto_anns$x,
                         Exp1 = log2(as.numeric(hto_sel[i,])+1),
                         Exp2 = log2(as.numeric(hto_sel[j,])+1))
    
    # log scaled
    gp1 <- ggplot(exp_df, mapping = aes(x = Exp1, y = Exp2, color = Annotation)) +
      geom_point(alpha = 0.5, size = 0.5) + xlab(paste("Log2 Read Counts", rownames(hto_sel)[i], sep = " ")) +
      ylab(paste("Log2 Read Counts", rownames(hto_sel)[j], sep = " "))
    plot(gp1)
    #dev.off()
  }
}
dev.off()
}
