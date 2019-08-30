library(Seurat)
library(Matrix)
library(DoubletDecon)
library(fields)
library(KernSmooth)
library(ROCR)
library(modes)

#modify the project name, e.g

#pr_names=c("CZI.PBMC","PBMC.8.HTO","Four.Cell.12.HTO")
pr_names=c("CZI.PBMC")

#for (j in 1:length(pr_names) ){

for (j in 1:length(pr_names) ){
    setwd(paste("/projects/ucar-lab/danaco/bncmrk-dblts/DoubletDecon/input/",pr_names[j],sep=""))
   
    #fl <- list.files(pattern = "preprocessed.seurat.for.DF.output.Rds")
    pr_name=pr_names[j]
       
    
    if (pr_names[j]=="CZI.PBMC"){
        #for (i in 1:10){
        for (i in c(6,8)){
            smpl=as.character(i)
            
location=paste0(getwd(),"/")
expressionFile=paste0(location,paste(pr_name,smpl,"preprocessed.seurat.for.DoubDec.counts.txt",sep="."))
expressionFile=gsub(x=expressionFile,pattern="\\.\\.",replacement="\\.")
genesFile=paste0(location, paste(pr_name,smpl,"preprocessed.seurat.for.DoubDec.Top50Genes.txt",sep="."))
genesFile=gsub(x=genesFile,pattern="\\.\\.",replacement="\\.")
clustersFile=paste0(location, paste(pr_name,smpl,"preprocessed.seurat.for.DoubDec.Cluster.txt",sep="."))
clustersFile=gsub(x=clustersFile,pattern="\\.\\.",replacement="\\.")        
write(expressionFile,"expressionFiles.txt",append=TRUE)
write(genesFile,"genesFiles.txt",append=TRUE)
write(clustersFile,"clustersFiles.txt",append=TRUE)}
}
    else {smpl=""
          
    location=paste0(getwd(),"/")
expressionFile=paste0(location,paste(pr_name,smpl,"preprocessed.seurat.for.DoubDec.counts.txt",sep="."))
expressionFile=gsub(x=expressionFile,pattern="\\.\\.",replacement="\\.")
genesFile=paste0(location, paste(pr_name,smpl,"preprocessed.seurat.for.DoubDec.Top50Genes.txt",sep="."))
genesFile=gsub(x=genesFile,pattern="\\.\\.",replacement="\\.")
clustersFile=paste0(location, paste(pr_name,smpl,"preprocessed.seurat.for.DoubDec.Cluster.txt",sep="."))
clustersFile=gsub(x=clustersFile,pattern="\\.\\.",replacement="\\.")
         write(expressionFile,"expressionFiles.txt",append=TRUE)
write(genesFile,"genesFiles.txt",append=TRUE)
write(clustersFile,"clustersFiles.txt",append=TRUE)}
}