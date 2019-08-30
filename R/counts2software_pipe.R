
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
                              numSamples=1){
  
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

scr_name2=paste(scr_name,"/output",sep="")
df_name2=paste(df_name,"/output",sep="")
doubdec_name2=paste(doubdec_name,"/output",sep="")
if (!foldersPresent){
system(paste("mkdir",scr_name,sep=" "))
system(paste("mkdir",df_name,sep=" "))
system(paste("mkdir",doubdec_name,sep=" "))
    
system(paste("mkdir",scr_name2,sep=" "))
system(paste("mkdir",df_name2,sep=" "))
system(paste("mkdir",doubdec_name2,sep=" "))
}
scr_name=paste(scr_name,"/",pr_name,sep="")
df_name=paste(df_name,"/",pr_name,sep="")
doubdec_name=paste(doubdec_name,"/",pr_name,sep="")
    
scr_name2=paste(scr_name2,"/",pr_name,sep="")
df_name2=paste(df_name2,"/",pr_name,sep="")
doubdec_name2=paste(doubdec_name2,"/",pr_name,sep="")
if (!foldersPresent){
system(paste("mkdir",scr_name,sep=" "))
system(paste("mkdir",df_name,sep=" "))
system(paste("mkdir",doubdec_name,sep=" "))
    
system(paste("mkdir",scr_name2,sep=" "))
system(paste("mkdir",df_name2,sep=" "))
system(paste("mkdir",doubdec_name2,sep=" "))
}

    
    
#scr_name=paste(home_dir,"/Scrublet/input/",pr_name,sep="")
#df_name=paste(home_dir,"/DF/input/",pr_name,sep="")
#doubdec_name=paste(home_dir,"/DoubletDecon/input/",pr_name,sep="")

#adt_name=paste(home_dir,"/Data/",pr_name,"/",pr_name,".ADT.raw.pooled.combined.matrix.Rds",sep="")
        
rna_name=paste(home_dir,"/Data/",pr_name,"/",pr_name,".RNA.raw.matrix.sparse.Rds",sep="")
#e.g, CZI.PBMC.raw.RNA.matrix.sparse.Rds stored in home/Data/CZI.Fourth.run folder

#hto_name=paste(home_dir,"/Data/",pr_name,"/",pr_name,".HTO.raw.pooled.sample.classifications.GMM.csv",sep="")
print(rna_name)
# rna counts matrix
comb.counts <- readRDS(rna_name)

print(dim(comb.counts))
scr_name=paste(scr_name,"/",pr_name,sep="")
df_name=paste(df_name,"/",pr_name,sep="")
doubdec_name=paste(doubdec_name,"/",pr_name,sep="")

# write matrix to output for Scrublet
if (numSamples==1){
writeMM(comb.counts, paste(scr_name, "raw.counts.for.Scrublet.mtx", sep = "."))
pbmc=CreateSeuratObject(counts=comb.counts,assay="RNA")
saveRDS(pbmc,paste(df_name,"raw.seurat.Input.for.DF.Rds",sep="."))
saveRDS(pbmc,paste(doubdec_name,"raw.seurat.Input.for.DoubDec.Rds",sep=".") )
}
if (numSamples>1){
  for(i in 1:numSamples){
      if (i<10){
   assign(paste("cc", i, sep = ""), 
          comb.counts[,substr(colnames(comb.counts)[],1,2)==paste(i,"-",sep="")])
          }
      else {
   assign(paste("cc", i, sep = ""), 
          comb.counts[,substr(colnames(comb.counts)[],1,2)==as.character(i)])
   #assign(paste("hto", i, sep = ""), hto_sel[substr(rownames(hto_sel)[],1,2)==as.character(i),])
          }
   saveRDS(CreateSeuratObject(counts=eval(parse(text=paste("cc", i, sep = ""))),
                          assay="RNA"),
          paste(df_name,i,"raw.seurat.Input","for.DF.Rds",sep=".") )
   
   saveRDS(CreateSeuratObject(counts=eval(parse(text=paste("cc", i, sep = ""))),
                             assay="RNA"),
          paste(doubdec_name,i,"raw.seurat.Input","for.DoubDec.Rds",sep=".") )
    
   writeMM(eval(parse(text=paste("cc", i, sep = ""))), 
           paste(scr_name, i, "raw.counts", "for.Scrublet.mtx", sep = "."))       
   }
 }

}
