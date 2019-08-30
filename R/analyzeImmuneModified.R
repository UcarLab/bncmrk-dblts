rm(list = ls())
set.seed(100)
library(Rtsne)
library(ggplot2)
library(Seurat)
library(Matrix)
setwd("~/Desktop/Codes/alphaModifiedPlots/groundTruths/")
dat1_raw=read.table("./1-HTO.all.cell.barcodes.for.demuxlet.txt.demux.filtered.VCF.best",header=TRUE,fill=TRUE)
dat2_raw=read.table("./2-HTO.all.cell.barcodes.for.demuxlet.txt.demux.filtered.VCF.best.best",header=TRUE,fill=TRUE)
dat3_raw=read.table("./3-HTO.all.cell.barcodes.for.demuxlet.txt.demux.filtered.VCF.best",header=TRUE,fill=TRUE)
dat4_raw=read.table("./4-HTO.all.cell.barcodes.for.demuxlet.txt.demux.filtered.VCF.best",header=TRUE,fill=TRUE)
dat5_raw=read.table("./5-HTO.all.cell.barcodes.for.demuxlet.txt.demux.filtered.VCF.best",header=TRUE,fill=TRUE)
dat6_raw=read.table("./6-HT",header=TRUE,fill=TRUE)
dat7_raw=read.table("~/Downloads/7-RNA.demux.filtered.VCF.best",header=TRUE,fill=TRUE)
dat8_raw=read.table("~/Downloads/8-RNA.demux.filtered.VCF.best",header=TRUE,fill=TRUE)
dat9_raw=read.table("~/Downloads/9-RNA.demux.filtered.VCF.best",header=TRUE,fill=TRUE)
dat10_raw=read.table("~/Downloads/10-RNA.demux.filtered.VCF.best",header=TRUE,fill=TRUE)

immune_raw=rbind(dat1_raw,dat2_raw,dat3_raw,dat4_raw,dat5_raw,
                 dat6_raw,dat7_raw,dat8_raw,dat9_raw,dat10_raw)

SNG.BEST.LLK=immune_raw$SNG.NEXT.GUESS

SNG.BEST.LLK=immune_raw$SNG.NEXT.GUESS

immune=data.frame(NUM.READS=immune_raw$NUM.READS,DIFF.LLK.SNG.DBL=SNG.BEST.LLK
                  -immune_raw$DBL.BEST.LLK, DROPLET.TYPE=immune_raw$DROPLET.TYPE)


#let's take the singlets only and fit a model
immuneSNG=subset(immune, DROPLET.TYPE=="SNG")
sng.m=lm(DIFF.LLK.SNG.DBL~NUM.READS,data=immuneSNG)

#let's take the doublets only and fit a model
immuneDBL=subset(immune, DROPLET.TYPE=="DBL")
dbl.m=lm(DIFF.LLK.SNG.DBL~NUM.READS,data=immuneDBL)

#let's take the doublets only and fit a model
immuneAMB=subset(immune, DROPLET.TYPE=="AMB")
amb.m=lm(DIFF.LLK.SNG.DBL~NUM.READS,data=immuneAMB)



ggplot(aes(NUM.READS,DIFF.LLK.SNG.DBL,color=DROPLET.TYPE),data=immune)+
  geom_abline(slope=sng.m$coefficients[2],intercept=sng.m$coefficients[1])+
  geom_abline(slope=dbl.m$coefficients[2],intercept=dbl.m$coefficients[1])+
  geom_point()+ggtitle('Freemuxlet Decision Boundary for immune')
ggsave("immuneDecisionBoundary.png")

#let's estimate the NUM.READ densities for both singlets & doublets
#and let's take their singlet or doublet ratios
#i.e, estimated mixing proportions alpha.hat, 1-alpha.hat into account

alpha.sng.hat=dim(immuneSNG)[1]/dim(immune)[1]
alpha.dbl.hat=dim(immuneDBL)[1]/dim(immune)[1]
alpha.amb.hat=dim(immuneAMB)[1]/dim(immune)[1]

immune.dens=density(immune$NUM.READS)
sng.dens=density(immuneSNG$NUM.READS)
dbl.dens=density(immuneDBL$NUM.READS)
amb.dens=density(immuneAMB$NUM.READS)

png(filename="immuneDensityFull.png")
plot(sng.dens$x,sng.dens$y*alpha.sng.hat,col="blue",
     type="l",xlab="NUM.READS",ylab="Density",xlim=c(0,4000),ylim=c(0,1.7e-3),
     main="immune Data Weighted densities of Singlets and Doublets")
lines(dbl.dens$x,dbl.dens$y*alpha.dbl.hat,col="red")
lines(amb.dens$x,dbl.dens$y*alpha.amb.hat,col="green")
lines(immune.dens,pch=4,cex=2,col="black")
legend(2700,1e-3,c("Singlet","Doublet","Ambiguous","All"), lwd=c(5,2),
       col=c("blue","red","green","black"))
dev.off()

immune.1500=subset(immune, NUM.READS<1500)
immune.1500.SNG=subset(immune.1500, DROPLET.TYPE=="SNG")
immune.1500.DBL=subset(immune.1500, DROPLET.TYPE=="DBL")
immune.1500.AMB=subset(immune.1500, DROPLET.TYPE=="AMB")


ggplot(aes(NUM.READS,DIFF.LLK.SNG.DBL,color=DROPLET.TYPE),data=immune.1500)+
  geom_point()+ggtitle('Freemuxlet Decision Boundary for immune, Num.Reads<1500')
ggsave("immuneDecisionBoundary1500.png")

alpha.1500.sng.hat=dim(immune.1500.SNG)[1]/dim(immune.1500)[1]
alpha.1500.dbl.hat=dim(immune.1500.DBL)[1]/dim(immune.1500)[1]
alpha.1500.amb.hat=dim(immune.1500.AMB)[1]/dim(immune.1500)[1]

immune.1500.dens=density(immune.1500$NUM.READS)
sng.1500.dens=density(immune.1500.SNG$NUM.READS)
dbl.1500.dens=density(immune.1500.DBL$NUM.READS)
amb.1500.dens=density(immune.1500.DBL$NUM.READS)

png(filename="immuneDensity1500.png",width=600)
plot(sng.1500.dens$x,sng.1500.dens$y*alpha.1500.sng.hat,col="blue",
     type="l",xlab="NUM.READS",ylab="Density",ylim=c(0,2e-3),
     main="immune Data Weighted densities of Singlets and Doublets NUM.READS<1500")
lines(dbl.1500.dens$x,dbl.1500.dens$y*alpha.1500.dbl.hat,col="red")
lines(amb.1500.dens$x,amb.1500.dens$y*alpha.1500.amb.hat,col="green")
lines(immune.1500.dens,pch=4,cex=2)
legend(1e3,3e-4,c("Singlet","Doublet","Ambiguous","All"), lwd=c(5,2),
       col=c("blue","red","green","black"))
dev.off()

source("~/Desktop/Codes/ToyModel/HTO.GMM.analysis.R")

immune.true=subset(s_anns,Sample_Short!="None")
