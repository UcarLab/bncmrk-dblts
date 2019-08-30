
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
library(dat)
library(tidyverse)

Scrublet_2.5_dbl_scores=read.csv("doublet_scores_2.5.csv")$X0
comb.counts=readRDS("/projects/ucar-lab/danaco/Ground/comb.counts.RDS")
groundTruth=read.csv("groundTruth.csv")[,2:3]
groundTruth$BARCODE=as.character(groundTruth$BARCODE)
groundTruth$S_type=as.character(groundTruth$S_type)
#gr_ids=which(as.character(groundTruth$BARCODE) %in% colnames(comb.counts))
#Scrublet_2.5_dbl_scores=Scrublet_2.5_dbl_scores[gr_ids]
DemuxletDefault=read.csv("DemuxletDefaultLabels.csv")[,2:3]
DemuxletDefault$DefDemux_labels=as.character(DemuxletDefault$DROPLET.TYPE)
DemuxletDefault$DROPLET.TYPE=NULL
#DemuxletDefault %>% dplyr::rename(DefDemux_labels = DROPLET.TYPE)
DemuxletModified=read.csv("DemuxletModifiedLabels.csv")[,2:3]
DemuxletModified$ModDemux_labels=as.character(DemuxletModified$DROPLET.TYPE)
DemuxletModified$DROPLET.TYPE=NULL
df=read.csv("../Software/DF/Ground/doubletFinderGround.csv")[,2:3]
#DemuxletDefault %>% dplyr::rename(ModDemux_labels = DROPLET.TYPE)

#Scrublet_2.5_dbl_scores=read.csv("doublet_scores_2.5.csv")$X0
scr_2.5=data.frame(BARCODE=colnames(comb.counts),scr_2.5_scores=Scrublet_2.5_dbl_scores)
groundTruth=merge(groundTruth,scr_2.5,by.x="BARCODE",by.y="BARCODE")
groundTruth=merge(groundTruth,DemuxletDefault,by.x="BARCODE",by.y="BARCODE")
groundTruth=merge(groundTruth,DemuxletModified,by.x="BARCODE",by.y="BARCODE")
groundTruth=merge(groundTruth,df,by.x="BARCODE",by.y="BARCODE")

groundTruth$scr_2.5_labels=groundTruth$scr_2.5_scores
groundTruth$scr_2.5_labels[groundTruth$scr_2.5_scores>0.25]="DBL"
groundTruth$scr_2.5_labels[groundTruth$scr_2.5_scores<=0.25]="SNG"

gr.SNG=subset(groundTruth,S_type=="SNG")
gr.DBL=subset(groundTruth,S_type=="DBL")

gr_2.5_SNG=subset(groundTruth,scr_2.5_labels=="SNG")
gr_2.5_DBL=subset(groundTruth,scr_2.5_labels=="DBL")
gr_Dem_SNG=subset(groundTruth,DefDemux_labels=="SNG")
gr_Dem_DBL=subset(groundTruth,DefDemux_labels=="DBL")
gr_df_SNG=subset(groundTruth,Df_label=="SNG")
gr_df_DBL=subset(groundTruth,Df_label=="DBL")

int_cells <- Reduce(intersect, list(gr.DBL$BARCODE, gr_2.5_SNG$BARCODE, gr_Dem_DBL$BARCODE,gr_df_SNG$BARCODE ))
length(int_cells)
#sum(gr.SNG$scr_2.5_labels=="SNG")/dim(gr.SNG)[1]
#sum(as.character(groundTruth$S_type)==as.character(DemuxletModified$DROPLET.TYPE))/length(groundTruth$S_type)

sum(groundTruth$S_type==as.character(groundTruth$Df_label))/length(groundTruth$S_type)
sum(gr.DBL$Df_label=="DBL")#/dim(gr.DBL)[1]
dim(gr_df_DBL)[1]/dim(groundTruth)[1]

Scrublet_2.5_dbl_scores=read.csv("doublet_scores_2.5.csv")$X0

#Scrublet_2.5_dbl_scores=Scrublet_2.5_dbl_scores[gr_ids]


demuxSNG=subset(DemuxletDefault, DROPLET.TYPE=="SNG")

groundTruth=read.csv("groundTruth.csv")


gr_ids=which(groundTruth$BARCODE %in% colnames(comb.counts))

#write.csv(groundTruth$BARCODE,"groundTruthBarcodes.csv")









dim(groundTruth)


