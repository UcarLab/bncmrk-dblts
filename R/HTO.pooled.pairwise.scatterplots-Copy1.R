# Plot HTO barcodes in pairwise manner, raw counts
rm(list = ls())
library(Matrix)
library(Biobase)
name <- "CZI_HTO_pooled"
setwd("/Users/lawlon/Documents/CZI/CITE-seq/PBMC_2018/Pooled_Runs/HTO_Analysis/")
# read in all hto
hto <- readRDS("CZI.PMBC.HTO.pooled.HTO.matrix.Rds")
bars <- c("Bcell-TGATGGCCTATTGGG", "Mono-AGTAAGTTCAGCGTA", "Tcell-TTCCGCCTCTCTTTG", "control-GTCAACTCTTTAGCG")
hto_sel <- hto[bars, ]

# pairwise scatterplots
pdf(file = paste(name, "raw_counts_pairwise_scatterplots.pdf", sep = "_"), onefile = T)
for (i in 1:nrow(hto_sel)) {
  j = i + 1
  for (j in j:nrow(hto_sel)) {
    # raw
    plot(x = as.numeric(hto_sel[i,]), y = as.numeric(hto_sel[j,]),
         main = "HTO Raw Read Counts", xlab = rownames(hto_sel)[i], ylab = rownames(hto_sel)[j], cex = 0.75)
    
    # log scaled
    plot(x = log2(as.numeric(hto_sel[i,])+1), y = log2(as.numeric(hto_sel[j,])+1),
         main = "HTO Log2 Raw Read Counts", xlab = rownames(hto_sel)[i], ylab = rownames(hto_sel)[j], cex = 0.75)
    
  }
}
dev.off()

# load in hto annotations
hto_anns <- read.csv("CZI.PMBC.HTO.pooled.sample.classifications.GMM.csv",
                     header = T, check.names = F, stringsAsFactors = F, row.names = 2)
table(rownames(hto_anns) == colnames(hto_sel))


for (i in 1:nrow(hto_sel)) {
  j = i + 1
  for (j in j:nrow(hto_sel)) {
    jpeg(file = paste(name, i, "vs", j, "pairwise_scatterplots_ggplots.jpeg", sep = "_"))
    # create a dataframe for ggplot
    exp_df <- data.frame(Cell_Barcode = colnames(hto_sel),
                         Annotation = hto_anns$Sample,
                         Exp1 = log2(as.numeric(hto_sel[i,])+1),
                         Exp2 = log2(as.numeric(hto_sel[j,])+1))
    
    # log scaled
    gp1 <- ggplot(exp_df, mapping = aes(x = Exp1, y = Exp2, color = Annotation)) +
      geom_point(alpha = 0.5, size = 0.5) + xlab(paste("Log2 Read Counts", rownames(hto_sel)[i], sep = " ")) +
      ylab(paste("Log2 Read Counts", rownames(hto_sel)[j], sep = " "))
    plot(gp1)
    dev.off()
  }
}
