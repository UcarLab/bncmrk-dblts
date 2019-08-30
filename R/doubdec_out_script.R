library(Seurat)
library(Matrix)
library(DoubletDecon)
library(fields)
library(KernSmooth)
library(ROCR)
library(modes)

#Running from /home_dir
source("/projects/ucar-lab/danaco/bncmrk-dblts/R/doubdec_outpipe.R")
    
#modify the project name, e.g
#setwd("/projects/ucar-lab/danaco/bncmrk-dblts/DoubletDecon/input")

#pr_names=c("CZI.PBMC","PBMC.8.HTO","Four.Cell.12.HTO")

#pr_names=c("CZI.PBMC","Four.Cell.12.HTO")
pr_names=c("CZI.PBMC")

for (j in 1:length(pr_names) ){
    setwd(paste("/projects/ucar-lab/danaco/bncmrk-dblts/DoubletDecon/input/",pr_names[j],sep=""))
   
   
    if (pr_names[j]=="CZI.PBMC"){
        #for (i in 9:10){
        for (i in c(8)){
            doubdec_outpipe(pr_name=pr_names[j],smpl=as.character(i))
        }
    }
    else {doubdec_outpipe(pr_name=pr_names[j])}
}
