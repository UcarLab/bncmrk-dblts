# Read in demuxlet object

#setwd("/projects/ucar-lab/danaco/Ground")

#name <- "CZI.Fourth.Run.RNA.ADT.UMAP.Ground.DF.optimized.10per.expect"

# list of Seurat objects
#fl <- list.files(pattern = "seurat.output.Rds")
# extract file name
df_outcleanpipe=function(file_name){

df.meta=read.csv(file_name)

if (length(grep(file_name,pattern="/"))>0){
    home_dir=paste(strsplit(file_name,split="/")[[1]][1:(length(strsplit(file_name,split="/")[[1]])-4)],collapse="/")
    pr_name=unlist(lapply(strsplit(file_name,split="/"),"[[",length(strsplit(file_name,split="/")[[1]])-1))
#setwd(paste("/",home_dir,"/DF/input/",pr_name,sep=""))
    file_name=unlist(lapply(strsplit(file_name,split="/"),"[[",length(strsplit(file_name,split="/")[[1]])))
    file_dir=paste(home_dir,"/DF/output/",pr_name,sep="")
}
else {home_dir=paste(strsplit(getwd(),split="/")[[1]][1:(length(strsplit(getwd(),split="/")[[1]])-3)],collapse="/")
      pr_name=unlist(lapply(strsplit(getwd(),split="/"),"[[",length(strsplit(getwd(),split="/")[[1]])))
      file_dir=paste(home_dir,"/DF/output/",pr_name,sep="")}

file_name=strsplit(x=file_name,split=".csv")[[1]][1]

name <- unlist(lapply(strsplit(x = file_name, split = ".DF.metadata"),`[[`, 1))

name=paste(name,"barcodes.labels.csv",sep=".")

name=paste0(file_dir,"/",name)

write.csv( data.frame(BARCODE=df.meta$X,LABELS=df.meta[,8]), name )
 }


