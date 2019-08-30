library(Seurat)
library(Matrix)
#library(DoubletDecon)
library(fields)
library(KernSmooth)
library(ROCR)
library(modes)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(dat)
library(tidyverse)
library(VennDiagram)
library(plyr)


train_test_validate_split=function(df,tr_pc=.7,tst_pc=.15,vl_pc=.15){

colnames(df)=gsub(x=colnames(df),pattern="X",replacement="")
colnames(df)=gsub(x=colnames(df),pattern="\\.",replacement="\\-")

spec = c(train = tr_pc, test = tst_pc, validate = vl_pc)


df_name=deparse(substitute(df))

g = sample(cut(
  seq(ncol(df)), 
  ncol(df)*cumsum(c(0,spec)),
  labels = names(spec)
))

res = list(df[,g=="train"],df[,g=="validate"],df[,g=="test"])
res
#assign(paste(df_name,"train",sep="."),df[,g=="train"])
#assign(paste(df_name,"train",sep="."),GetAssayData(NormalizeData(CreateSeuratObject(eval(parse(text=paste(df_name,"train",sep=".")))),method="CLR")))
#assign(paste(df_name,"validate",sep="."),df[,g=="validate"])
#assign(paste(df_name,"validate",sep="."),GetAssayData(NormalizeData(CreateSeuratObject(eval(parse(text=paste(df_name,"validate",sep=".")))),method="CLR")))
#assign(paste(df_name,"test",sep="."),df[,g=="test"])
#assign(paste(df_name,"test",sep="."),GetAssayData(NormalizeData(CreateSeuratObject(eval(parse(text=paste(df_name,"test",sep=".")))),method="CLR")))
}



setwd("/projects/ucar-lab/danaco/bncmrk-dblts/DF/input/CZI.PBMC/")
fl=list.files(pattern="preprocessed.seurat.for.DF.output")
print("read fl")
for (i in 4:length(fl)){
    print(i)
    name=strsplit(fl[i],split=".preprocessed")[[1]][1]
    print(name)
    #print(strsplit(fl[i],split=".preprocessed")[[1]][1])
    res=train_test_validate_split(data.frame(GetAssayData(readRDS(fl[i]))))
    assign(paste0(name,".train"),res[1])
    write.csv(eval(substitute(paste0(name,".train"))),paste0(name,".train.csv"))

    assign(paste0(name,".validate"),res[2])
    write.csv(eval(substitute(paste0(name,".validate"))),paste0(name,".validate.csv"))

    assign(paste0(name,".test"),res[3])
    write.csv(eval(substitute(paste0(name,".test"))),paste0(name,".test.csv"))

}
