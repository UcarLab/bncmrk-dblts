library(Matrix)
library(Seurat)

HTO.MULTIseqDemux.analysis=function(hto_mtx){
  
  sr_hto <- CreateSeuratObject(counts = hto_mtx,assay="HTO")#, min.features = 300)
  sr_hto <- NormalizeData(sr_hto, assay = "HTO", normalization.method = "CLR")
  sr_hto <- MULTIseqDemux(sr_hto, assay = "HTO")
  sr_hto$MULTI_ID=as.character(sr_hto$MULTI_ID)
  sr_hto$MULTI_ID[sr_hto$MULTI_ID!="Negative"&sr_hto$MULTI_ID!="Doublet"]="Singlet"
 # multiseq_ids=unlist(lapply(sr_hto$MULTI_classification,str_count,pattern="_"))
  #sr_hto$MULTI_classification[multiseq_ids==0 & sr_hto$MULTI_classification!="Negative"]="Singlet"
  #sr_hto$MULTI_classification[multiseq_ids>0]="Doublet"
  #sr_hto$MULTI_classification[sr_hto$MULTI_classification!="Doublet"&sr_hto$MULTI_classification!="Negative"]="Singlet"
  #sr_hto$MULTI_classification=as.factor(sr_hto$MULTI_classification)
  sr_hto$MULTI_ID=as.factor(sr_hto$MULTI_ID)
  
  sr_hto

}