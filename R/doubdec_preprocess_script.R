library(Seurat)
library(Matrix)
library(DoubletFinder)
library(fields)
library(KernSmooth)
library(ROCR)
library(modes)

#Running from /home_dir
source("/projects/ucar-lab/danaco/bncmrk-dblts/R/doubdec_preprocess.R")


#modify the project name, e.g
pr_names=c("CZI.PBMC","PBMC.8.HTO","Four.Cell.12.HTO")

#pr_names=c("CZI.PBMC","Four.Cell.12.HTO")

for (j in 1:length(pr_names) ){
  setwd(paste("/projects/ucar-lab/danaco/bncmrk-dblts/DoubletDecon/input/",pr_names[j],sep=""))
  fl <- list.files(pattern = "raw.seurat.Input.for.DoubDec.Rds")
  for (i in 1:length(fl)){
    doubdec_preprocess(fl[i])}
  }