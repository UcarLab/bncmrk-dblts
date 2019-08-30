# Read in demuxlet object
library(Seurat)
library(Matrix)
library(DoubletFinder)
library(fields)
library(KernSmooth)
library(ROCR)
library(modes)
#setwd("/projects/ucar-lab/danaco/Ground")

#name <- "CZI.Fourth.Run.RNA.ADT.UMAP.Ground.DF.optimized.10per.expect"

# list of Seurat objects
#fl <- list.files(pattern = "seurat.output.Rds")
# extract file name
df_outpipe=function(file_name){

seu_kidney=readRDS(file_name)


if (length(grep(file_name,pattern="/"))>0){
    home_dir=paste(strsplit(file_name,split="/")[[1]][1:(length(strsplit(file_name,split="/")[[1]])-4)],collapse="/")
    pr_name=unlist(lapply(strsplit(file_name,split="/"),"[[",length(strsplit(file_name,split="/")[[1]])-1))
#setwd(paste("/",home_dir,"/DF/input/",pr_name,sep=""))
    file_name=unlist(lapply(strsplit(file_name,split="/"),"[[",length(strsplit(file_name,split="/")[[1]])))
    file_dir=paste(home_dir,"/DF/output/",pr_name,sep="")
}
else {home_dir=paste(strsplit(getwd(),split="/")[[1]][1:(length(strsplit(getwd(),split="/")[[1]])-3)],collapse="/")
      pr_name=unlist(lapply(strsplit(getwd(),split="/"),"[[",length(strsplit(getwd(),split="/")[[1]])))
      file_dir=paste(home_dir,"/DF/output/",pr_name,sep="")}

file_name=strsplit(x=file_name,split=".Rds")[[1]][1]

name <- unlist(lapply(strsplit(x = file_name, split = ".preprocessed"),`[[`, 1))

name=paste(name,"DF.output",sep=".")

# read through list and	do analyses

 
print(paste("Computation commences for",name))

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
print("Parameter Sweep in progress...")
sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
#saveRDS(sweep.res.list_kidney,"try_sweep.res.list_kidney.Rds")
print("Summarize Sweep in progress...")

sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
#saveRDS(sweep.stats_kidney,"try_sweep.res.list_kidney.Rds")
print("Peak finding in progress...")

bcmvn_kidney <- find.pK(sweep.stats_kidney)
#saveRDS(bcmvn_kidney,"try_bcmvn_kidney.Rds")

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
#sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
#gt.calls <- seu_kidney@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"]
#sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
#bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
print("Extracting annotations...")
annotations <- seu_kidney@meta.data$ClusteringResults
#saveRDS(annotations,"try_annotations.Rds")
print("Finding homotypic proportions...")
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
print("Estimating exp_poi...")
nExp_poi <- round(0.075*length(as.character(seu_kidney@meta.data$RNA_snn_res.1.2)))  ## Assuming 7.5% doublet formation rate - #tailor for your dataset
print("Rounding nExp_poi and homotypic...")
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#write.csv(paste(homotypic.prop,nExp_poi,nExp_poi.adj,sep=", "),"try_homotypic.csv") 
    
 ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
 # identify is the maxima BCM value
print("Identifying the maxima BCM Value...")
 id_max <- which(bcmvn_kidney$BCmetric == max(bcmvn_kidney$BCmetric))
 # extract the pK value
print("Extracting the pK value...")
 pk_use <- as.numeric(as.character(bcmvn_kidney$pK[id_max]))    
#write.csv(paste(id_max,pk_use,sep=", "),"try_maximaBCMpK.RData")


## Run DoubletFinder with varying classification stringencies -----------------------
print("Running doublet finder without re-used pANNs")
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = pk_use, nExp = nExp_poi, reuse.pANN = FALSE,)
#saveRDS(seu_kidney,"try_seu_kidney_without_pANN.Rds")
  #  seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = #FALSE)


 #ser <- doubletFinder_v3(ser, PCs = 1:20, pN = 0.25, pK = pk_use, nExp = nExp_poi, reuse.pANN = FALSE)
 # extract metadata
print("Extracting meta data")
 meta_data <- seu_kidney@meta.data
#saveRDS(meta_data,"try_meta.data1.Rds")
 # determine which column is pANN
print("Determining which column is pANN")
 id_pann <- which(grepl(x = colnames(meta_data), pattern = "pANN"))
#write.csv(id_pann,"try_id_paNN.csv")
    
print("Running doublet finder again by re-using known pANNs")    
  seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = pk_use, nExp = nExp_poi.adj, reuse.pANN = colnames(meta_data)[id_pann])
# saveRDS(seu_kidney,"try_seu_kidney_with_pANN.Rds")
 # identify first classifications
print("Identifying first classifications. First extracting metadata...")
 meta_data <- seu_kidney@meta.data
#saveRDS(meta_data,"try_meta_data2.Rds")
print("Identifying first classifications. Second extracting DF.classifications from metadata...")

 id_cl <- which(grepl(x = colnames(meta_data), pattern = "DF.classifications"))
 #write.csv(id_cl,"try_id_cl.csv")
    
#seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "#pANN_0.25_0.09_913", sct = FALSE)
 
print("Adjusting hi_lo DF...")    
df_labels <- seu_kidney@meta.data[,id_cl[1]]
df_barcodes= colnames(GetAssayData(seu_kidney))
 #print("Adjusting doublets, singlets, doublet_lo s...")
 #seu_kidney@meta.data$DF_hi.lo[which(seu_kidney@meta.data$DF_hi.lo == "Doublet" & seu_kidney@meta.data[, id_cl[2]] == "Singlet")] <- "Doublet_lo"
 #printing("Adjusting doublets, doublet hi s...")
 #seu_kidney@meta.data$DF_hi.lo[which(seu_kidney@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
    
# save new seurat object
 print("Saving the processed seurat object to Rds...")
 saveRDS(seu_kidney, file = paste(paste(file_dir,name,sep="/"), "doubletfinder.seurat", pk_use, "Rds", sep = "."))
 
# make a table of doubletfinder classifications
print("Saving the table of DF classifications...")
 write.csv(table(df_labels, df_barcodes), file = paste(paste(file_dir,name,sep="/"), "doublet.finder.numbers", pk_use, "csv", sep = "."))

 # table of cells and their id
print("Saving the table of cells and their id s...")
 write.csv(seu_kidney@meta.data, file = paste(paste(file_dir,name,sep="/"), "DF.metadata", pk_use, "csv", sep = "."))
 
    
}


