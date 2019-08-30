    
set.seed(100)
library(Rtsne)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(Matrix)

df_preprocess_seurat=function(file_name){

pbmc= readRDS(file_name)
        
#if the pipeline is run from the /home_dir/DF/input/pr_name folder
#with file_name="Foo.Rds",then that's ok. don't need to change anything
#but run from another folder with locations
#e.g, file_name="/home_dir/DF/input/pr_name/Foo.Rds", then dissect it
#to get home folder, file_name, pr_name etc such that
#you get a meaningful output on /home_dir/DF/input/pr_name folder
    
if (length(grep(file_name,pattern="/"))>0){
    home_dir=paste(strsplit(file_name,split="/")[[1]][1:(length(strsplit(file_name,split="/")[[1]])-4)],collapse="/")
    pr_name=unlist(lapply(strsplit(file_name,split="/"),"[[",length(strsplit(file_name,split="/")[[1]])-1))
#setwd(paste("/",home_dir,"/DF/input/",pr_name,sep=""))
    file_name=unlist(lapply(strsplit(file_name,split="/"),"[[",length(strsplit(file_name,split="/")[[1]])))
    file_dir=paste(home_dir,"/DF/input/",pr_name,sep="")
}
else {file_dir=getwd()}
    
#setwd(file_dir)   
#Take the raw *raw.seurat.Input.for.DF.Rds file, which are preprocessed seurat objects
#after passing it through this pipeline
#normalize->find features->scale->pca->umap->neighbors->clusters
#the processed output files are ready to go to next doublet finder steps
    
    
name=unlist(lapply(strsplit(file_name,split=".raw"),"[[",1))

name=paste(name,"preprocessed.seurat.for.DF",sep=".")



#pbmc=CreateSeuratObject(counts=comb.counts,meta.data=hto_sel,assay="RNA")
#print("Creating seurat object")
    
pdf(file = paste(paste(file_dir,name,sep="/"), "plots.pdf", sep = "."), onefile = T)


## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
pbmc <- NormalizeData(object = pbmc, normalization.method = "CLR")
pbmc <- FindVariableFeatures(object = pbmc, selection.method = 'vst', nfeatures = 500)
pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc))
print("Running PCA")
pbmc= RunPCA(pbmc,features = VariableFeatures(object = pbmc), verbose = FALSE)
print("Running UMAP")
pbmc <- RunUMAP(object = pbmc, dims = 1:10)
pbmc <- FindNeighbors(object = pbmc, dims = 1:10)
pbmc <- FindClusters(object = pbmc, resolution = 1.2)




# save seurat object to R data file
saveRDS(pbmc, file = paste(paste(file_dir,name,sep="/"), "output.Rds", sep = "."))
# scale the data
dev.off()
}
