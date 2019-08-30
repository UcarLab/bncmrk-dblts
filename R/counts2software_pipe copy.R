
set.seed(100)
library(Rtsne)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(Matrix)

counts2software_pipe=function(pr_name, home_dir=getwd(),foldersPresent=TRUE,
                              geneIds=FALSE, adtPresent=FALSE, numSamples=1,
                               rw_sm=50,cl_sm=400){
  
#pr_name : use projects name such as CZI.PBMC 
#if you have DF, Scrublet, DoubletDecon (with sub)folders, modify foldersPresent=TRUE
#if you are not running this function in the home folder, modify home_dir

#our pbmc-hto project has 10 samples, modify it to 10
#our gene ids are given, turn that option to TRUE for our data
#our data has ADTs, turn adtPresent to TRUE
 
#Default for RNA matrix cleaning, modify if needed
#If a gene's expressed less than 50 times within all the droplets, rowSums < 50 eliminate gene from RNA count
#If a droplet contains less than 400 genes, colSums <400 eliminate droplet from RNA count
   
    
# change directory to base/home (grandest ancestry) directory
# R codes are stored on /home_dir/R
# Data on home_dir/Data/pr_name etc etc

#if (grep(file_name,pattern="/"){home_dir=unlist(lapply(strsplit(file_name,split="/"),"[[",2))}
   
 
setwd(home_dir)
    
scr_name=paste(home_dir,"/Scrublet",sep="")
df_name=paste(home_dir,"/DF",sep="")
doubdec_name=paste(home_dir,"/DoubletDecon",sep="")
if (!foldersPresent){
system(paste("mkdir",scr_name,sep=" "))
system(paste("mkdir",df_name,sep=" "))
system(paste("mkdir",doubdec_name,sep=" "))
}
scr_name=paste(scr_name,"/input",sep="")
df_name=paste(df_name,"/input",sep="")
doubdec_name=paste(doubdec_name,"/input",sep="")
if (!foldersPresent){
system(paste("mkdir",scr_name,sep=" "))
system(paste("mkdir",df_name,sep=" "))
system(paste("mkdir",doubdec_name,sep=" "))
}
scr_name=paste(scr_name,"/",pr_name,sep="")
df_name=paste(df_name,"/",pr_name,sep="")
doubdec_name=paste(doubdec_name,"/",pr_name,sep="")
if (!foldersPresent){
system(paste("mkdir",scr_name,sep=" "))
system(paste("mkdir",df_name,sep=" "))
system(paste("mkdir",doubdec_name,sep=" "))
}

    
    
#scr_name=paste(home_dir,"/Scrublet/input/",pr_name,sep="")
#df_name=paste(home_dir,"/DF/input/",pr_name,sep="")
#doubdec_name=paste(home_dir,"/DoubletDecon/input/",pr_name,sep="")

#adt_name=paste(home_dir,"/Data/",pr_name,"/",pr_name,".ADT.raw.pooled.combined.matrix.Rds",sep="")
        
rna_name=paste(home_dir,"/Data/",pr_name,"/",pr_name,".RNA.raw.matrix.sparse.Rds",sep="")
#e.g, CZI.PBMC.raw.RNA.matrix.sparse.Rds stored in home/Data/CZI.Fourth.run folder
        
hto_name=paste(home_dir,"/Data/",pr_name,"/",pr_name,".HTO.raw.pooled.sample.classifications.GMM.csv",sep="")

               

                
# rna counts matrix
all_rna <- readRDS(rna_name)
# remove genes with zeros across all cells (non-informative)
all_rna <- all_rna[rowSums(all_rna) > 0, ]

if (geneIds){

genes_name=paste(home_dir,"/Data/",pr_name,"/",pr_name,".genes.tsv",sep="")

# read in genes
genes <- read.delim(genes_name,
                    header = F, stringsAsFactors = F, check.names = F)
# change rownames of genes from ENSEMBL id to symbols
g_ids <- NULL
for (g in 1:nrow(all_rna)) {
  idx <- which(genes$V1 == rownames(all_rna)[g])
  g_ids <- c(g_ids, genes$V2[idx[1]])
}
rownames(all_rna) <- g_ids

}

# remove duplicate gene symbols from rna matrix
dup_id <- which(duplicated(rownames(all_rna)))
all_rna <- all_rna[-dup_id, ]

#Clear out the genes not expressed enough, and droplets not containing enough genes
all_rna=all_rna[rowSums(all_rna)>rw_sm,colSums(all_rna)>cl_sm]

# read file containing HTO identity info (singlet, doublet, etc)
hto <- read.csv(hto_name,
                header = T, check.names = F, stringsAsFactors = F, row.names = 2)
#change column names of HTO matrix, the prefixes of the columns in the HTO matrix need to match those of the RNA matrix
rownames(hto) <- gsub(x = rownames(hto), pattern = "HTO", replacement = "Sample")
# should have ~300k plus cell barcodes
length(which(colnames(all_rna) %in% rownames(hto)))

if (adtPresent){
    
adt_name=paste(home_dir,"/Data/",pr_name,"/",pr_name,".raw.ADT.pooled.combined.matrix.Rds",sep="")
    
# read in adt/protein matrix data, this will also be used to cluster data
adt_counts <- readRDS(adt_name)
colnames(adt_counts) <- gsub(x = colnames(adt_counts), pattern = "ADT", replacement = "Sample")

# some rows in the ADT matrix contain summary info and not levels of specific proteins, so remove them
adt_id <- which(rownames(adt_counts) %in% c("bad_struct", "no_match", "total_reads"))
adt_counts <- adt_counts[-adt_id,]
# remove control ids too (these aren't useful)
adt_cont <- which(grepl(x = rownames(adt_counts), pattern = "control"))
adt_counts <- adt_counts[-adt_cont, ]
rownames(adt_counts) <- paste("CITE", rownames(adt_counts), sep = "-")

# change underscores in ADT row names to hyphens
rownames(adt_counts) <- gsub(x = rownames(adt_counts), pattern = "_", replacement = "-")
    
}
# analysis for all cells
hto_sel <- hto

#Taking out the empty droplets from HTO None s
no_none_barcodes=rownames(hto[as.character(hto$Sample) != "None",]) #cells that doesn't have none barcode
cells_w_nb=Reduce(intersect,list(colnames(all_rna),no_none_barcodes)) #cells with no none barcodes overlapping with rna
no_none_ids=which(colnames(all_rna) %in% cells_w_nb) #where do these cells located in rna matrix
all_rna=all_rna[,no_none_ids] #Taking those cells out of all_rna matrix

if (adtPresent){
# find intersection of all cell barcodes across HTO, RNA, and ADT matrices
int_cells <- Reduce(intersect, list(rownames(hto_sel), colnames(all_rna), colnames(adt_counts)))
rna_ids <- which(colnames(all_rna) %in% int_cells)
rna_sel <- all_rna[, rna_ids]
adt_sel <- adt_counts[, colnames(rna_sel)]
hto_sel <- hto[colnames(rna_sel), ]
    
# rename metadata column names
#colnames(hto_sel) <- c("ID", "Markers", "Treatment", "Run")

write.csv(rownames(adt_sel),paste(paste(home_dir,"Data",pr_name,"/")
                                  ,"adt.features.csv",sep="."))
    
# combine rna-seq and adt into a single matrix
comb.counts <- Matrix::rbind2(rna_sel, as.matrix(adt_sel))
}
else {
# find intersection of all cell barcodes across HTO, RNA, and ADT matrices
int_cells <- Reduce(intersect, list(rownames(hto_sel), colnames(all_rna)))
rna_ids <- which(colnames(all_rna) %in% int_cells)
rna_sel <- all_rna[, rna_ids]
hto_sel <- hto[colnames(rna_sel), ]
    
# rename metadata column names
#colnames(hto_sel) <- c("ID", "Markers", "Treatment", "Run")
    
# combine rna-seq and adt into a single matrix
comb.counts <- rna_sel
}
# do some filtering of the gene matrix (genes with low counts aren't useful for clustering)
#rna_sel <- rna_sel[Matrix::rowSums(rna_sel) > 5,]
#rna_sel[rowSums(rna_sel)>40,colSums(rna_sel)>400]

scr_name=paste(scr_name,"/",pr_name,sep="")
df_name=paste(df_name,"/",pr_name,sep="")
doubdec_name=paste(doubdec_name,"/",pr_name,sep="")

# write matrix to output for Scrublet
writeMM(comb.counts, paste(scr_name, "raw.counts.for.Scrublet.mtx", sep = "."))
pbmc=CreateSeuratObject(counts=comb.counts,meta.data=hto_sel,assay="RNA")
saveRDS(pbmc,paste(df_name,"raw.seurat.Input.for.DF.Rds",sep="."))
saveRDS(pbmc,paste(doubdec_name,"raw.seurat.Input.for.DoubDec.Rds",sep=".") )
 if (numSamples>1){
  for(i in 1:numSamples){
      if (ii<10){
   assign(paste("cc", i, sep = ""), comb.counts[,substr(colnames(comb.counts)[],1,2)==paste(i,"-",sep="")]
   assign(paste("hto", i, sep = ""), hto_sel[substr(rownames(hto_sel)[],1,2)==paste(i,"-",sep=""),])
          }
      else {
   assign(paste("cc", i, sep = ""), comb.counts[,substr(colnames(comb.counts)[],1,2)==as.character(i)]
   assign(paste("hto", i, sep = ""), hto_sel[substr(rownames(hto_sel)[],1,2)==as.character(i),])
          }
   saveRDS(CreateSeuratObject(counts=eval(parse(text=quote(paste("cc", i, sep = "")))),
                             meta.data=eval(parse(text=quote(paste("hto", i, sep = "")))),
                             assay="RNA"),
          paste(df_name,i,"raw.seurat.Input","for.DF.Rds",sep=".") )
   
   saveRDS(CreateSeuratObject(counts=eval(parse(text=quote(paste("cc", i, sep = "")))),
                             meta.data=eval(parse(text=quote(paste("hto", i, sep = "")))),
                             assay="RNA"),
          paste(doubdec_name,i,"raw.seurat.Input","for.DoubDec.Rds",sep=".") )
    
   writeMM(eval(parse(text=quote(paste("cc", i, sep = "")))), 
           paste(scr_name, i, "raw.counts", "for.Scrublet.mtx", sep = "."))       
   }
 }

}
