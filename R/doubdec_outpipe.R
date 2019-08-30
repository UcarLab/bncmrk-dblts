# Read in demuxlet object
library(Seurat)
library(Matrix)
library(DoubletDecon)
library(fields)
library(KernSmooth)
library(ROCR)
library(modes)
#setwd("/projects/ucar-lab/danaco/Ground")

#name <- "CZI.Fourth.Run.RNA.ADT.UMAP.Ground.DF.optimized.10per.expect"

# list of Seurat objects
#fl <- list.files(pattern = "seurat.output.Rds")
# extract file name
doubdec_outpipe=function(pr_name,smpl=""){

source("/projects/ucar-lab/danaco/bncmrk-dblts/R/Seurat_Pre_Process2.R")
#ser=readRDS(file_name)


#if (length(grep(file_name,pattern="/"))>0){
#    home_dir=paste(strsplit(file_name,split="/")[[1]][1:(length(strsplit(file_name,split="/")[[1]])-4)],collapse="/")
#    pr_name=unlist(lapply(strsplit(file_name,split="/"),"[[",length(strsplit(deli,split="/")[[1]])-1))
#setwd(paste("/",home_dir,"/DF/input/",pr_name,sep=""))
#    file_name=unlist(lapply(strsplit(file_name,split="/"),"[[",length(strsplit(deli,split="/")[[1]])))
#    file_dir=paste(home_dir,"/DoubletDecon/out/",pr_name,sep="")
#}
#else {
      home_dir=paste(strsplit(getwd(),split="/")[[1]][1:(length(strsplit(getwd(),split="/")[[1]])-3)],collapse="/")
      #pr_name=unlist(lapply(strsplit(getwd(),split="/"),"[[",length(strsplit(getwd(),split="/")[[1]])))
      file_dir=paste(home_dir,"/DoubletDecon/output/",pr_name,sep="")

#file_name=strsplit(x=file_name,split=".Rds")[[1]][1]

#name <- unlist(lapply(strsplit(x = file_name, split = ".preprocessed"),`[[`, 1))

#name=paste(name,"DF.output",sep=".")

# read through list and	do analyses

#name=paste(pr_name

location=paste0(getwd(),"/")
print( paste("Reading from location",location) )
expressionFile=paste0(location,paste(pr_name,smpl,"preprocessed.seurat.for.DoubDec.counts.txt",sep="."))
expressionFile=gsub(x=expressionFile,pattern="\\.\\.",replacement="\\.")
print( paste("Reading expressionFile",expressionFile) )
genesFile=paste0(location, paste(pr_name,smpl,"preprocessed.seurat.for.DoubDec.Top50Genes.txt",sep="."))
genesFile=gsub(x=genesFile,pattern="\\.\\.",replacement="\\.")
print( paste("Reading genesFile",genesFile) )
clustersFile=paste0(location, paste(pr_name,smpl,"preprocessed.seurat.for.DoubDec.Cluster.txt",sep="."))
clustersFile=gsub(x=clustersFile,pattern="\\.\\.",replacement="\\.")
print( paste("Reading clustersFile",clustersFile) )


newFiles=Seurat_Pre_Process2(expressionFile, genesFile, clustersFile)
    
filename=paste(pr_name,smpl,sep=".")
location=paste(file_dir,"",sep="/")
    
write.table(newFiles$newExpressionFile, paste0(location, filename, "_expression"), sep="\t")
write.table(newFiles$newFullExpressionFile, paste0(location, filename, "_fullExpression"), sep="\t")
write.table(newFiles$newGroupsFile, paste0(location, filename , "_groups"), sep="\t", col.names = F)

pdf(paste0(location,filename,"_heatmaps.pdf"))
results=Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                           groupsFile=newFiles$newGroupsFile, 
                           filename=filename, 
                           location=location,
                           fullDataFile=NULL, 
                           removeCC=FALSE, 
                           species="hsa", 
                           rhop=2.75, 
                           write=TRUE, 
                           PMF=TRUE, 
                           useFull=FALSE, 
                           heatmap=FALSE,
                           centroids=TRUE,
                           num_doubs=100, 
                           only50=FALSE,
                           min_uniq=4)    
dev.off()
}


