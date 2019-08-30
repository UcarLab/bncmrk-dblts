library(Matrix)
library(Seurat)

HTO.HTODemux.analysis=function(hto_mtx){

sr_hto <- CreateSeuratObject(counts = hto_mtx,assay="HTO")#, min.features = 300)
sr_hto <- NormalizeData(sr_hto, assay = "HTO", normalization.method = "CLR")
sr_hto <- HTODemux(sr_hto, assay = "HTO", positive.quantile = 0.99)
sr_hto
}