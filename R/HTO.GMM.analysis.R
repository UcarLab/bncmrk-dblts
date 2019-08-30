# load in CITE-seq data
library(data.table)
library(Matrix)
library(Seurat)
library(mclust)
#rm(list = ls())

HTO.GMM.analysis=function(hto_mtx){
# CLR normalize the data for GMM modeling
sr_obj <- CreateSeuratObject(counts =hto_mtx )
mrkr_num=dim(hto_mtx)[1]
print(mrkr_num)
# centered log-normalization
sr_obj <- NormalizeData(sr_obj, normalization.method = "CLR") # norm data in @data slot
#norm_counts <- sr_obj@data
norm_counts <- GetAssayData(sr_obj, slot = "data")
#rownames(norm_counts) <- c("Control", "Mono", "Bcell", "Tcell")
                             
# use gmm method to classify cell types
p.exp <- norm_counts
p.bin <- p.exp
p.bin[] <- 0
p.names <- rownames(p.exp)
                             
# keep list of output paramters
params <- list()
for(i in 1:nrow(p.exp))
{
  print(i)
 x <- Mclust(p.exp[i,],G=2)
p.bin[i,] <- x$classification
params[[i]] <- x$parameters
names(params)[i] <- rownames(norm_counts)[i]
  }
                             
                             
                             c.sum <- apply(p.bin,2,sum)
                             sample.type <- character(length = ncol(p.bin))
                             s.type <- character(length = ncol(p.bin))
                             tr.s.type <- character(length = ncol(p.bin))
                             
                             print(length(sample.type[c.sum < mrkr_num]))
                             # mrkr_num markers, if sum is mrkr_num, then sample cannot be classified
                             sample.type[c.sum == mrkr_num] <- "Negative"
                             s.type[c.sum == mrkr_num] <- "Negative"
                             tr.s.type[c.sum == mrkr_num] <- "Negative"
                             
                             for(i in 1:length(c.sum)) {
                               if(c.sum[i] ==mrkr_num+1 & length(which(p.bin[,i] ==2) == 1)) {
                                 sample.type[i] <- p.names[which(p.bin[,i] ==2)]
                                 s.type[i]="Singlet"
                                 tr.s.type[i]="Singlet"
                               }
                               
                               if(c.sum[i] > mrkr_num+1 & length(which(p.bin[,i] ==2) == 2)) {
                                 ids.sel <- p.names[which(p.bin[,i] ==2)]
                                 nams <- ""
                                 for (num in 1:length(ids.sel)) {
                                   nams <- paste(ids.sel[num], nams, sep="+")  
                                 }
                                 nams <- gsub(x = nams, pattern = "\\++$", replacement = "")
                                 sample.type[i] <- nams
                               }
                               
                               if(c.sum[i] > mrkr_num+1 & length(which(p.bin[,i] ==2)) == 2) {
                                 s.type[i]="Doublet"
                                 tr.s.type[i]="Doublet"
                               } else if (c.sum[i] >mrkr_num+1  & length(which(p.bin[,i] ==2)) == 3){
                                 s.type[i]="Doublet"
                                 tr.s.type[i]="TRP"
                               } else if (c.sum[i] > mrkr_num+1 & length(which(p.bin[,i] ==2)) == 4){
                                 s.type[i]="Doublet" #hto12.htos=HTO.GMM.analysis(hto12.htos)
                                 tr.s.type[i]="QDR"
                               } else if (c.sum[i] > mrkr_num+1 & length(which(p.bin[,i] ==2)) == 5){
                                 s.type[i]="Doublet" #hto12.htos=HTO.GMM.analysis(hto12.htos)
                                 tr.s.type[i]="QNT"
                               }else if (c.sum[i] > mrkr_num+1 & length(which(p.bin[,i] ==2)) == 6){
                                 s.type[i]="Doublet" #hto12.htos=HTO.GMM.analysis(hto12.htos)
                                 tr.s.type[i]="SXT"
                               }else if (c.sum[i] > mrkr_num+1 & length(which(p.bin[,i] ==2)) == 7){
                                 s.type[i]="Doublet" #hto12.htos=HTO.GMM.analysis(hto12.htos)
                                 tr.s.type[i]="SPT"
                               }
                               
                             }
                             
                             # short sample assignment
                             d_ids <- which(grepl(x = sample.type, pattern = "\\+"))
                             Sample_Short <- sample.type
                             Sample_Short[d_ids] <- "Doublet"
                             
                             # data frame of cell barcodes and assignments
                             s_anns <- data.frame(Cell_Barcode = colnames(p.bin), Sample = sample.type, 
                                                  Sample_Short = Sample_Short, S_type=s.type
                                                  ,TR_S_type=tr.s.type)
                             table(s_anns$Sample)
                             table(s_anns$Sample_Short)
                             # add run identifiers
                             runs <- strsplit(x = as.character(s_anns$Cell_Barcode), split = "_")
                             run_nam <- unlist(lapply(runs, `[[`, 1))
                             s_anns$Run <- run_nam
                             # proportions of each sample type
                             tab_s <- table(s_anns$Run, s_anns$Sample_Short)
                             rsums <- rowSums(tab_s)
                             per_s <- 100 * (tab_s / rsums)
                             
                          s_anns
            }