set.seed(100)
library(Rtsne)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(Matrix)


doubdec_preprocess=function(file_name){
    
source("/projects/ucar-lab/danaco/bncmrk-dblts/R/Top50Markers.R")
    
#if (length(grep(file_name,pattern="/"))>0){setwd(unlist(lapply(strsplit(file_name,split="/"),"[[",2)))}

pbmc= readRDS(file_name)

print(paste("Processing",file_name))
    
    
#if the pipeline is run from the /home_dir/DF/input/pr_name folder
#with file_name="Foo.Rds",then set the directory to /home_dir/DoubletDecon/input/pr_name
#but if run from another folder with locations given with "/"
#e.g, file_name="/home_dir/DF/input/pr_name/Foo.Rds", then dissect it
#to get home folder, file_name, pr_name etc such that
#you get a meaningful output on /home_dir/DoubletDecon/input/pr_name folder
    
if (length(grep(file_name,pattern="/"))>0){
    home_dir=paste(strsplit(file_name,split="/")[[1]][1:(length(strsplit(file_name,split="/")[[1]])-4)],collapse="/")
    pr_name=unlist(lapply(strsplit(file_name,split="/"),"[[",length(strsplit(file_name,split="/")[[1]])-1))
#setwd(paste("/",home_dir,"/DoubletDecon/input/",pr_name,sep=""))
    file_name=unlist(lapply(strsplit(file_name,split="/"),"[[",length(strsplit(file_name,split="/")[[1]])))
    file_dir=paste(home_dir,"/DoubletDecon/input/",pr_name,sep="")

}
else {home_dir=paste(strsplit(getwd(),split="/")[[1]][1:(length(strsplit(getwd(),split="/")[[1]])-3)],collapse="/")
      pr_name=unlist(lapply(strsplit(getwd(),split="/"),"[[",length(strsplit(getwd(),split="/")[[1]])))
      file_dir=paste(home_dir,"/DoubletDecon/input/",pr_name,sep="")}

#setwd(file_dir)
    
name=unlist(lapply(strsplit(file_name,split=".raw"),"[[",1))

name=paste(paste(file_dir,name,sep="/"),"preprocessed.seurat.for.DoubDec",sep=".")

pbmc <- NormalizeData(object = pbmc, normalization.method = "CLR")
pbmc <- FindVariableFeatures(object = pbmc, selection.method = 'vst', nfeatures = 500)
pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc))
print("Running PCA")
pbmc= RunPCA(pbmc,features = VariableFeatures(object = pbmc), verbose = FALSE)
print("Running UMAP")
pbmc <- RunUMAP(object = pbmc, dims = 1:10)
pbmc <- FindNeighbors(object = pbmc, dims = 1:10)
pbmc <- FindClusters(object = pbmc, resolution = 1.2)


    
print(paste("Writing them to name", name))

write.table(GetAssayData(pbmc),paste(name,"counts.txt",sep="."),sep="\t", col.names=NA)
    
write.table(pbmc@active.ident,paste(name,"Cluster.txt",sep="."),sep="\t",col.names=NA)    

cm=FindAllMarkers(pbmc) #cluster marker info
    
write.csv(cm,paste(name,"AllMarkers.csv",sep="."))

topmrkrs=Top50Markers(cm)

write.table(topmrkrs,paste(name,"Top50Genes.txt",sep="."),sep="\t") 

#write.table(top_n(n=50,x=mrkrs,wt="avg_logFC"),paste(name,"Top50Genes.txt",sep="."),sep="\t")    
    
    
}
