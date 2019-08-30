
setwd("~/bncmrk-dblts")
source("~/bncmrk-dblts/R/counts2software_pipe.R")

counts2software_pipe(pr_name="CZI.PBMC",foldersPresent=FALSE,numSamples = 10)
counts2software_pipe(pr_name="PBMC.8.HTO",foldersPresent=FALSE,numSamples = 1)
counts2software_pipe(pr_name="Four.Cell.12.HTO",foldersPresent=FALSE,numSamples = 1)
