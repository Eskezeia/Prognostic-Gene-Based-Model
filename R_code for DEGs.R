# Identification of Differential expressed genes(DEGs)
library(limma)
library(edgeR)
setwd("wdpath") # working directroy 
dat1<-read.csv("inputfileName",row.names=1) # Normalized gene expression matrix{ Refer TCGA dataset the link given below}
dim(dat1)
dat2<-log2(dat1 + 1)

dat3<-as.matrix(dat2)
dim(dat3)
design.mat<-read.csv("design_HCC_gene.csv") # matrix samples containing tissue type (Normal and Tumor)
design.mat<-design.mat[,-1]
design.mat
dim(design.mat)
contrast.mat<-matrix(c(1,-1), ncol = 1)
dimnames(contrast.mat)<-list(c('T', 'N'), "Diff")
sample<-factor(rep(c("T", "N"), c(371,50)))
design.mat<-model.matrix(~0+sample)
colnames(design.mat)<-levels(sample)
design.mat
dim(design.mat)
contrast.mat<-makeContrasts(diff = T - N, levels = design.mat)
contrast.mat
fit<-lmFit(dat3, design.mat)
fit2<-contrasts.fit(fit, contrast.mat)
fit3<-eBayes(fit2)
deg<-topTable(fit3) # top ranked DEGs
























