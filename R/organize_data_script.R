library(Seurat)
library(Matrix)

#First read PBMC.8.HTO data
# Load in the UMI matrix

pbmc.umis=readRDS("~/bncmrk-dblts/Data/PBMC.8.HTO/pbmc_umi_mtx.rds")

# For generating a hashtag count matrix from FASTQ files, please refer to
# https://github.com/Hoohm/CITE-seq-Count.  Load in the HTO count matrix
pbmc.htos <- readRDS("~/bncmrk-dblts/Data/PBMC.8.HTO/pbmc_hto_mtx.rds")

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
#rownames(pbmc.htos)

# Read Four.cells data

hto12.htos <- readRDS("~/bncmrk-dblts/Data/Four.Cell.12.HTO/hto12_hto_mtx.rds")
hto12.umis <- readRDS("~/bncmrk-dblts/Data/Four.Cell.12.HTO/hto12_umi_mtx.rds")

cells.use <- intersect(rownames(hto12.htos), colnames(hto12.umis))

hto12.htos=t(hto12.htos[cells.use,1:12])
hto12.umis=hto12.umis[,cells.use]


#Read CZI.PBMC data

all_rna=readRDS("~/bncmrk-dblts/Data/CZI.PBMC/CZI.PBMC.RNA.raw.matrix.sparse.unprocessed.Rds")

genes=read.delim("~/bncmrk-dblts/Data/CZI.PBMC/CZI.PBMC.genes.tsv",
                          header = F, stringsAsFactors = F, check.names = F)

# change rownames of genes from ENSEMBL id to symbols
g_ids <- NULL
for (g in 1:nrow(all_rna)) {
  idx <- which(genes$V1 == rownames(all_rna)[g])
  g_ids <- c(g_ids, genes$V2[idx[1]])
}
rownames(all_rna) <- g_ids

all_rna=all_rna[order(rownames(all_rna)),]
# remove duplicate gene symbols from rna matrix
dup_id <- which(duplicated(rownames(all_rna)))
all_rna <- all_rna[-dup_id, ]

all_hto=readRDS("~/bncmrk-dblts/Data/CZI.PBMC/CZI.PBMC.HTO.ordered.count.matrix.Rds")

#all_hto[c("control-GTCAACTCTTTAGCG", "Mono-AGTAAGTTCAGCGTA", 
         # "Bcell-TGATGGCCTATTGGG", "Tcell-TTCCGCCTCTCTTTG"),]

all_hto=all_hto[c(1,2,3,5),]
#rownames(all_hto) <- c("Control", "Mono", "Bcell", "Tcell")
rownames(all_hto) <- c("Bcell", "Mono", "Tcell", "Control")

#Clean all the data sets
pbmc.umis=pbmc.umis[rowSums(pbmc.umis)>0,colSums(pbmc.umis)>0]
hto12.umis=hto12.umis[rowSums(hto12.umis)>0,colSums(hto12.umis)>0]
all_rna=all_rna[rowSums(all_rna)>0,colSums(all_rna)>0]

all_hto=all_hto[,intersect(colnames(all_hto),colnames(all_rna))]
source("~/bncmrk-dblts/R/HTO.GMM.analysis.R")
source("~/bncmrk-dblts/R/HTO.HTODemux.analysis.R")
source("~/bncmrk-dblts/R/HTO.MULTIseqDemux.analysis.R")

pbmc.htos.gmm=HTO.GMM.analysis(pbmc.htos)
write.csv(pbmc.htos.gmm[c(1,4)]
,"~/bncmrk-dblts/GroundTruths/PBMC.8.HTO/PBMC.8.HTO.HTO.sample.classifications.GMM.csv")
pbmc.htos.htod=HTO.HTODemux.analysis(pbmc.htos)
write.csv(pbmc.htos.htod$HTO_classification.global,
"~/bncmrk-dblts/GroundTruths/PBMC.8.HTO/PBMC.8.HTO.HTO.sample.classifications.HTODemux.csv")
pbmc.htos.msq=HTO.MULTIseqDemux.analysis(pbmc.htos)
write.csv(pbmc.htos.msq,
          "~/bncmrk-dblts/GroundTruths/PBMC.8.HTO/PBMC.8.HTO.HTO.sample.classifications.MULTIseqDemux.csv")

hto12.htos.gmm=HTO.GMM.analysis(hto12.htos)
write.csv(hto12.htos.gmm[c(1,4)]
          ,"~/bncmrk-dblts/GroundTruths/Four.Cell.12.HTO/Four.Cell.12.HTO.HTO.sample.classifications.GMM.csv")
hto12.htos.htod=HTO.HTODemux.analysis(hto12.htos)
write.csv(hto12.htos.htod$HTO_classification.global,
          "~/bncmrk-dblts/GroundTruths/Four.Cell.12.HTO/Four.Cell.12.HTO.HTO.sample.classifications.HTODemux.csv")
hto12.htos.msq=HTO.MULTIseqDemux.analysis(hto12.htos)
write.csv(hto12.htos.msq,
          "~/bncmrk-dblts/GroundTruths/Four.Cell.12.HTO/Four.Cell.12.HTO.HTO.sample.classifications.MULTIseqDemux.csv")

library(Seurat)
library(Matrix)

#First read PBMC.8.HTO data
# Load in the UMI matrix

pbmc.umis=readRDS("~/bncmrk-dblts/Data/PBMC.8.HTO/pbmc_umi_mtx.rds")

# For generating a hashtag count matrix from FASTQ files, please refer to
# https://github.com/Hoohm/CITE-seq-Count.  Load in the HTO count matrix
pbmc.htos <- readRDS("~/bncmrk-dblts/Data/PBMC.8.HTO/pbmc_hto_mtx.rds")

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
#rownames(pbmc.htos)

# Read Four.cells data

hto12.htos <- readRDS("~/bncmrk-dblts/Data/Four.Cell.12.HTO/hto12_hto_mtx.rds")
hto12.umis <- readRDS("~/bncmrk-dblts/Data/Four.Cell.12.HTO/hto12_umi_mtx.rds")

cells.use <- intersect(rownames(hto12.htos), colnames(hto12.umis))

hto12.htos=t(hto12.htos[cells.use,1:12])
hto12.umis=hto12.umis[,cells.use]


#Read CZI.PBMC data

all_rna=readRDS("~/bncmrk-dblts/Data/CZI.PBMC/CZI.PBMC.RNA.raw.matrix.sparse.unprocessed.Rds")

genes=read.delim("~/bncmrk-dblts/Data/CZI.PBMC/CZI.PBMC.genes.tsv",
                 header = F, stringsAsFactors = F, check.names = F)

# change rownames of genes from ENSEMBL id to symbols
g_ids <- NULL
for (g in 1:nrow(all_rna)) {
  idx <- which(genes$V1 == rownames(all_rna)[g])
  g_ids <- c(g_ids, genes$V2[idx[1]])
}
rownames(all_rna) <- g_ids

all_rna=all_rna[order(rownames(all_rna)),]
# remove duplicate gene symbols from rna matrix
dup_id <- which(duplicated(rownames(all_rna)))
all_rna <- all_rna[-dup_id, ]

all_hto=readRDS("~/bncmrk-dblts/Data/CZI.PBMC/CZI.PBMC.HTO.ordered.count.matrix.Rds")

#all_hto[c("control-GTCAACTCTTTAGCG", "Mono-AGTAAGTTCAGCGTA", 
# "Bcell-TGATGGCCTATTGGG", "Tcell-TTCCGCCTCTCTTTG"),]

all_hto=all_hto[c(1,2,3,5),]
#rownames(all_hto) <- c("Control", "Mono", "Bcell", "Tcell")
rownames(all_hto) <- c("Bcell", "Mono", "Tcell", "Control")

#Clean all the data sets
pbmc.umis=pbmc.umis[rowSums(pbmc.umis)>0,colSums(pbmc.umis)>0]
hto12.umis=hto12.umis[rowSums(hto12.umis)>0,colSums(hto12.umis)>0]
all_rna=all_rna[rowSums(all_rna)>0,colSums(all_rna)>0]

all_hto=all_hto[,intersect(colnames(all_hto),colnames(all_rna))]
source("~/bncmrk-dblts/R/HTO.GMM.analysis.R")
source("~/bncmrk-dblts/R/HTO.HTODemux.analysis.R")
source("~/bncmrk-dblts/R/HTO.MULTIseqDemux.analysis.R")

pbmc.htos.gmm=HTO.GMM.analysis(pbmc.htos)
write.csv(pbmc.htos.gmm[c(1,4)]
          ,"~/bncmrk-dblts/GroundTruths/PBMC.8.HTO/PBMC.8.HTO.HTO.sample.classifications.GMM.csv")
write.csv(pbmc.htos.gmm[c(1,2)]
          ,"~/bncmrk-dblts/GroundTruths/PBMC.8.HTO/PBMC.8.HTO.HTO.sample.annotations.GMM.csv")
pbmc.htos.htod=HTO.HTODemux.analysis(pbmc.htos)
write.csv(pbmc.htos.htod$HTO_classification.global,
          "~/bncmrk-dblts/GroundTruths/PBMC.8.HTO/PBMC.8.HTO.HTO.sample.classifications.HTODemux.csv")
write.csv(pbmc.htos.htod$HTO_classification,
          "~/bncmrk-dblts/GroundTruths/PBMC.8.HTO/PBMC.8.HTO.HTO.sample.annotations.HTODemux.csv")
pbmc.htos.msq=HTO.MULTIseqDemux.analysis(pbmc.htos)
write.csv(pbmc.htos.msq$MULTI_ID,
          "~/bncmrk-dblts/GroundTruths/PBMC.8.HTO/PBMC.8.HTO.HTO.sample.classifications.MULTIseqDemux.csv")
write.csv(pbmc.htos.msq$MULTI_classification,
          "~/bncmrk-dblts/GroundTruths/PBMC.8.HTO/PBMC.8.HTO.HTO.sample.annotations.MULTIseqDemux.csv")

hto12.htos.gmm=HTO.GMM.analysis(hto12.htos)
write.csv(hto12.htos.gmm[c(1,4)]
          ,"~/bncmrk-dblts/GroundTruths/Four.Cell.12.HTO/Four.Cell.12.HTO.HTO.sample.classifications.GMM.csv")
write.csv(hto12.htos.gmm[c(1,2)]
          ,"~/bncmrk-dblts/GroundTruths/Four.Cell.12.HTO/Four.Cell.12.HTO.HTO.sample.annotations.GMM.csv")
hto12.htos.htod=HTO.HTODemux.analysis(hto12.htos)
write.csv(hto12.htos.htod$HTO_classification.global,
          "~/bncmrk-dblts/GroundTruths/Four.Cell.12.HTO/Four.Cell.12.HTO.HTO.sample.classifications.HTODemux.csv")
write.csv(hto12.htos.htod$HTO_classification,
          "~/bncmrk-dblts/GroundTruths/Four.Cell.12.HTO/Four.Cell.12.HTO.HTO.sample.annotations.HTODemux.csv")
hto12.htos.msq=HTO.MULTIseqDemux.analysis(hto12.htos)
write.csv(hto12.htos.msq$MULTI_ID,
          "~/bncmrk-dblts/GroundTruths/Four.Cell.12.HTO/Four.Cell.12.HTO.HTO.sample.classifications.MULTIseqDemux.csv")
write.csv(hto12.htos.msq$MULTI_classification,
          "~/bncmrk-dblts/GroundTruths/Four.Cell.12.HTO/Four.Cell.12.HTO.HTO.sample.annotations.MULTIseqDemux.csv")

czi.pbmc.gmm=HTO.GMM.analysis(all_hto)
write.csv(czi.pbmc.gmm[c(1,4)]
          ,"~/bncmrk-dblts/GroundTruths/CZI.PBMC/CZI.PBMC.HTO.sample.classifications.GMM.csv")
write.csv(czi.pbmc.gmm[c(1,2)]
          ,"~/bncmrk-dblts/GroundTruths/CZI.PBMC/CZI.PBMC.HTO.sample.annotations.GMM.csv")
czi.pbmc.htod=HTO.HTODemux.analysis(all_hto)
write.csv(czi.pbmc.htod$HTO_classification.global,
          "~/bncmrk-dblts/GroundTruths/CZI.PBMC/CZI.PBMC.HTO.sample.classifications.HTODemux.csv")
write.csv(czi.pbmc.htod$HTO_classification,
          "~/bncmrk-dblts/GroundTruths/CZI.PBMC/CZI.PBMC.HTO.sample.annotations.HTODemux.csv")
czi.pbmc.msq=HTO.MULTIseqDemux.analysis(all_hto)
write.csv(czi.pbmc.msq$MULTI_ID,
          "~/bncmrk-dblts/GroundTruths/CZI.PBMC/CZI.PBMC.HTO.sample.classifications.MULTIseqDemux.csv")
write.csv(czi.pbmc.msq$MULTI_classification,
          "~/bncmrk-dblts/GroundTruths/CZI.PBMC/CZI.PBMC.HTO.sample.annotations.MULTIseqDemux.csv")


#genes to be used in benchmarking and neural network projects
genes.nn=Reduce(intersect,list(rownames(hto12.umis),rownames(pbmc.umis),rownames(all_rna)))
write.csv(genes.nn,"/projects/ucar-lab/danaco/bncmrk-dblts/Data/genes.nn.csv")

pbmc.umis=pbmc.umis[genes.nn,]
hto12.umis=hto12.umis[genes.nn,]
all_rna=all_rna[genes.nn,]


saveRDS(all_rna,"/projects/ucar-lab/danaco/bncmrk-dblts/Data/CZI.PBMC/CZI.PBMC.RNA.raw.matrix.sparse.Rds")
saveRDS(pbmc.umis,"/projects/ucar-lab/danaco/bncmrk-dblts/Data/PBMC.8.HTO/PBMC.8.HTO.RNA.raw.matrix.sparse.Rds")
saveRDS(hto12.umis,"/projects/ucar-lab/danaco/bncmrk-dblts/Data/Four.Cell.12.HTO/Four.Cell.12.HTO.RNA.raw.matrix.sparse.Rds")

pbmc.gmm=read.csv("/projects/ucar-lab/danaco/bncmrk-dblts/GroundTruths/PBMC.8.HTO/PBMC.8.HTO.HTO.sample.classifications.GMM.csv")[,2:3]
pbmc.gmm_nnbrcds=pbmc.gmm$Cell_Barcode[pbmc.gmm$S_type!="Negative"]
write.csv(pbmc.gmm_nnbrcds,"/projects/ucar-lab/danaco/bncmrk-dblts/GroundBarcodes/PBMC.8.HTO/PBMC.8.HTO.non.neg.brcds.GMM.csv")
hto12.gmm=read.csv("/projects/ucar-lab/danaco/bncmrk-dblts/GroundTruths/Four.Cell.12.HTO/Four.Cell.12.HTO.HTO.sample.classifications.GMM.csv")[,2:3]
hto12.gmm_nnbrcds=hto12.gmm$Cell_Barcode[hto12.gmm$S_type!="Negative"]
write.csv(pbmc.gmm_nnbrcds,"/projects/ucar-lab/danaco/bncmrk-dblts/GroundBarcodes/Four.Cell.12.HTO/Four.Cell.12.HTO.non.neg.brcds.GMM.csv")
all_rna.gmm=read.csv("/projects/ucar-lab/danaco/bncmrk-dblts/GroundTruths/CZI.PBMC/CZI.PBMC.HTO.sample.classifications.GMM.csv")[,2:3]
all_rna.gmm_nnbrcds=all_rna.gmm$Cell_Barcode[all_rna.gmm$S_type!="Negative"]
write.csv(all_rna.gmm_nnbrcds,"/projects/ucar-lab/danaco/bncmrk-dblts/GroundBarcodes/CZI.PBMC/CZI.PBMC.non.neg.brcds.GMM.csv")

pbmc.htdm=read.csv("/projects/ucar-lab/danaco/bncmrk-dblts/GroundTruths/PBMC.8.HTO/PBMC.8.HTO.HTO.sample.classifications.HTODemux.csv")
pbmc.htdm_nnbrcds=pbmc.htdm$X[pbmc.htdm$x!="Negative"]
write.csv(pbmc.htdm_nnbrcds,"/projects/ucar-lab/danaco/bncmrk-dblts/GroundBarcodes/PBMC.8.HTO/PBMC.8.HTO.non.neg.brcds.HTODemux.csv")
hto12.htdm=read.csv("/projects/ucar-lab/danaco/bncmrk-dblts/GroundTruths/Four.Cell.12.HTO/Four.Cell.12.HTO.HTO.sample.classifications.HTODemux.csv")
hto12.htdm_nnbrcds=hto12.htdm$X[hto12.htdm$x!="Negative"]
write.csv(pbmc.htdm_nnbrcds,"/projects/ucar-lab/danaco/bncmrk-dblts/GroundBarcodes/Four.Cell.12.HTO/Four.Cell.12.HTO.non.neg.brcds.HTODemux.csv")
all_rna.htdm=read.csv("/projects/ucar-lab/danaco/bncmrk-dblts/GroundTruths/CZI.PBMC/CZI.PBMC.HTO.sample.classifications.HTODemux.csv")
all_rna.htdm_nnbrcds=all_rna.htdm$X[all_rna.htdm$x!="Negative"]
write.csv(all_rna.htdm_nnbrcds,"/projects/ucar-lab/danaco/bncmrk-dblts/GroundBarcodes/CZI.PBMC/CZI.PBMC.non.neg.brcds.HTODemux.csv")

pbmc.msdm=read.csv("/projects/ucar-lab/danaco/bncmrk-dblts/GroundTruths/PBMC.8.HTO/PBMC.8.HTO.HTO.sample.classifications.MULTIseqDemux.csv")
pbmc.msdm_nnbrcds=pbmc.msdm$X[pbmc.msdm$x!="Negative"]
write.csv(pbmc.msdm_nnbrcds,"/projects/ucar-lab/danaco/bncmrk-dblts/GroundBarcodes/PBMC.8.HTO/PBMC.8.HTO.non.neg.brcds.MULTIseqDemux.csv")
hto12.msdm=read.csv("/projects/ucar-lab/danaco/bncmrk-dblts/GroundTruths/Four.Cell.12.HTO/Four.Cell.12.HTO.HTO.sample.classifications.MULTIseqDemux.csv")
hto12.msdm_nnbrcds=hto12.msdm$X[hto12.msdm$x!="Negative"]
write.csv(pbmc.msdm_nnbrcds,"/projects/ucar-lab/danaco/bncmrk-dblts/GroundBarcodes/Four.Cell.12.HTO/Four.Cell.12.HTO.non.neg.brcds.MULTIseqDemux.csv")
all_rna.msdm=read.csv("/projects/ucar-lab/danaco/bncmrk-dblts/GroundTruths/CZI.PBMC/CZI.PBMC.HTO.sample.classifications.MULTIseqDemux.csv")
all_rna.msdm_nnbrcds=all_rna.msdm$X[all_rna.msdm$x!="Negative"]
write.csv(all_rna.msdm_nnbrcds,"/projects/ucar-lab/danaco/bncmrk-dblts/GroundBarcodes/CZI.PBMC/CZI.PBMC.non.neg.brcds.MULTIseqDemux.csv")



