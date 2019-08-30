    
set.seed(100)
library(Rtsne)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(Matrix)

df_preprocess_seurat=function(file_name,features_name){

pbmc= readRDS(file_name)
        
#if the pipeline is run from the /home_dir/DF/input/pr_name folder
#with file_name="Foo.Rds",then that's ok. don't need to change anything
#but run from another folder with locations
#e.g, file_name="/home_dir/DF/input/pr_name/Foo.Rds", then dissect it
#to get home folder, file_name, pr_name etc such that
#you get a meaningful output on /home_dir/DF/input/pr_name folder
    
if (grep(file_name,pattern="/"){
    home_dir=paste(strsplit(file_name,split="/")[[1]][1:(length(strsplit(file_name,split="/")[[1]])-4)],collapse="/")
    pr_name=unlist(lapply(strsplit(file_name,split="/"),"[[",length(strsplit(deli,split="/")[[1]])-1))
#setwd(paste("/",home_dir,"/DF/input/",pr_name,sep=""))
    file_name=unlist(lapply(strsplit(file_name,split="/"),"[[",length(strsplit(deli,split="/")[[1]])))
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

if(!missing(features_name)) {
ftrs_lst=read.csv(features_name)[,2]
}


#pbmc=CreateSeuratObject(counts=comb.counts,meta.data=hto_sel,assay="RNA")
#print("Creating seurat object")
    
pdf(file = paste(paste(file_dir,name,sep="/"), "plots.pdf", sep = "."), onefile = T)
pbmc <- NormalizeData(object = pbmc, normalization.method = "CLR")

# identify the top 500 variable features to be used for dimension reduction and clustering
pbmc <- FindVariableFeatures(object = pbmc, selection.method = 'vst', nfeatures = 500)


print(length(x = VariableFeatures(object = pbmc)))

# scale the data
pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc))
print("Normalizing and scaling")

# perform principal component analysis 
#pbmc <- RunPCA(object = pbmc, features = union(rownames(adt_sel), VariableFeatures(object = pbmc)), verbose = FALSE)

if(!missing(features_name)) {
pbmc <- RunPCA(object = pbmc, features = union(ftrs_lst, VariableFeatures(object = pbmc)), verbose = FALSE)
}

print("Running PCA")

# perform UMAP and identify clusters based on RNA matrix
pbmc <- RunUMAP(object = pbmc, dims = 1:20)
pbmc <- FindNeighbors(object = pbmc, dims = 1:20)
pbmc <- FindClusters(object = pbmc, resolution = 1.2)

# cell proportions by cluster
tab <- table(pbmc$RNA_snn_res.1.2, pbmc$Treatment)
print(tab)
write.csv(x = tab, file = paste(paste(file_dir,name,sep="/"), "cell.proportions.csv", sep = "."))

# output the metadata for each single cell
meta_d <- pbmc@meta.data
write.csv(x = meta_d, file = paste(paste(file_dir,name,sep="/"), "cell.metadata.csv", sep = "."))

# save seurat object to R data file
saveRDS(pbmc, file = paste(paste(file_dir,name,sep="/"), "output.Rds", sep = "."))
# scale the data
dev.off()
}