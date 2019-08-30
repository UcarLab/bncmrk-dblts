
rm(list = ls())
set.seed(100)
library(Rtsne)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(Matrix)

ftrs_lst=read.csv("../../../Ground/groundTruthBarcodes.csv")[,2]
df.1.md=read.csv("1_.DF.output.DF.metadata.0.3.csv")
pres_ids=which(as.character(df.1.md$X) %in% ftrs_lst)
df.1.md=df.1.md[pres_ids,]
df.10.md=read.csv("10_.DF.output.DF.metadata.0.01.csv")
pres_ids=which(as.character(df.10.md$X) %in% ftrs_lst)
df.10.md=df.10.md[pres_ids,]
df.2.md=read.csv("2_.DF.output.DF.metadata.0.01.csv")
pres_ids=which(as.character(df.2.md$X) %in% ftrs_lst)
df.2.md=df.2.md[pres_ids,]
df.3.md=read.csv("3_.DF.output.DF.metadata.0.16.csv")
pres_ids=which(as.character(df.3.md$X) %in% ftrs_lst)
df.3.md=df.3.md[pres_ids,]
df.4.md=read.csv("4_.DF.output.DF.metadata.0.3.csv")
pres_ids=which(as.character(df.4.md$X) %in% ftrs_lst)
df.4.md=df.4.md[pres_ids,]
df.5.md=read.csv("5_.DF.output.DF.metadata.0.18.csv")
pres_ids=which(as.character(df.5.md$X) %in% ftrs_lst)
df.5.md=df.5.md[pres_ids,]
df.6.md=read.csv("6_.DF.output.DF.metadata.0.25.csv")
pres_ids=which(as.character(df.6.md$X) %in% ftrs_lst)
df.6.md=df.6.md[pres_ids,]
df.7.md=read.csv("7_.DF.output.DF.metadata.0.26.csv")
pres_ids=which(as.character(df.7.md$X) %in% ftrs_lst)
df.7.md=df.7.md[pres_ids,]
df.8.md=read.csv("8_.DF.output.DF.metadata.0.27.csv")
pres_ids=which(as.character(df.8.md$X) %in% ftrs_lst)
df.8.md=df.8.md[pres_ids,]
df.9.md=read.csv("9_.DF.output.DF.metadata.0.02.csv")
pres_ids=which(as.character(df.9.md$X) %in% ftrs_lst)
df.9.md=df.9.md[pres_ids,]

dr_1=sum(as.character(df.1.md$DF.classifications_0.25_0.3_2325)=="Doublet")
dr_10=sum(as.character(df.10.md$DF.classifications_0.25_0.01_1913)=="Doublet")
dr_2=sum(as.character(df.2.md$DF.classifications_0.25_0.01_2269)=="Doublet")
dr_3=sum(as.character(df.3.md$DF.classifications_0.25_0.16_2164)=="Doublet")
dr_4=sum(as.character(df.4.md$DF.classifications_0.25_0.3_2208)=="Doublet")
dr_5=sum(as.character(df.5.md$DF.classifications_0.25_0.18_2139)=="Doublet")
dr_6=sum(as.character(df.6.md$DF.classifications_0.25_0.25_2135)=="Doublet")
dr_7=sum(as.character(df.7.md$DF.classifications_0.25_0.26_2206)=="Doublet")
dr_8=sum(as.character(df.8.md$DF.classifications_0.25_0.27_2056)=="Doublet")
dr_9=sum(as.character(df.9.md$DF.classifications_0.25_0.02_2274)=="Doublet")

df.1=data.frame(BARCODE=as.character(df.1.md$X),Df_label=as.character(df.1.md$DF.classifications_0.25_0.3_2325))
df.10=data.frame(BARCODE=as.character(df.10.md$X),Df_label=as.character(df.10.md$DF.classifications_0.25_0.01_1913))
df.2=data.frame(BARCODE=as.character(df.2.md$X),Df_label=as.character(df.2.md$DF.classifications_0.25_0.01_2269))
df.3=data.frame(BARCODE=as.character(df.3.md$X),Df_label=as.character(df.3.md$DF.classifications_0.25_0.16_2164))
df.4=data.frame(BARCODE=as.character(df.4.md$X),Df_label=as.character(df.4.md$DF.classifications_0.25_0.3_2208))
df.5=data.frame(BARCODE=as.character(df.5.md$X),Df_label=as.character(df.5.md$DF.classifications_0.25_0.18_2139))
df.6=data.frame(BARCODE=as.character(df.6.md$X),Df_label=as.character(df.6.md$DF.classifications_0.25_0.25_2135))
df.7=data.frame(BARCODE=as.character(df.7.md$X),Df_label=as.character(df.7.md$DF.classifications_0.25_0.26_2206))
df.8=data.frame(BARCODE=as.character(df.8.md$X),Df_label=as.character(df.8.md$DF.classifications_0.25_0.27_2056))
df.9=data.frame(BARCODE=as.character(df.9.md$X),Df_label=as.character(df.9.md$DF.classifications_0.25_0.02_2274))

df.gr=rbind(df.1,df.2,df.3,df.4,df.5,df.6,df.7,df.8,df.9,df.10)
df.gr$Df_label=as.character(df.gr$Df_label)

df.gr$Df_label[df.gr$Df_label=="Singlet"]="SNG"
df.gr$Df_label[df.gr$Df_label=="Doublet"]="DBL"

saveRDS(df.gr,"doubletFinderGround.Rds")
write.csv(df.gr,"doubletFinderGround.csv")

#df.gr$Df_label[which(as.character(df.gr$Df_label)=="Doublet")]






