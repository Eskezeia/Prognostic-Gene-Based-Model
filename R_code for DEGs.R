# I am running this code and it works. Save the data in one file and setwd in the same folder folder.
# do step by step, if there is any problem inform me or may be   we discuss via skype, Thenk you.
# save the data attached blelow and do DEG analysis.

library(limma)
library(edgeR)
setwd("C:/Users/Reon-Bruk/Desktop/geneT2/kidney2")
dat1<-read.csv("mRNA_KIRC.csv",row.names=1)
dim(dat1)
dat2<-log2(dat1 + 1)

dat3<-as.matrix(dat2)
dim(dat3)
design.mat<-read.csv("design_KIRC_gene.csv")
design.mat<-design.mat[,-1]
design.mat
dim(design.mat)
contrast.mat<-matrix(c(1,-1), ncol = 1)
dimnames(contrast.mat)<-list(c('T', 'N'), "Diff")
sample<-factor(rep(c("T", "N"), c(534,72)))
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
setwd("C:/Users/Reon-Bruk/Desktop/cmu_geneData")

dat1<-read.csv("HCC_FiveGene_TPM_survival2.csv")
write.csv(dat1, "dat.csv")
dat2<-read.csv("dat.csv")
head(dat2)

model<- coxph(Surv(os, event) ~ dat2$LCATx, data =dat2)
summary(model)
model<- coxph(Surv(y.time,y.status) ~ RSx, data =dat1)
summary(model)

fit2 <- survfit(Surv(y.time,y.status) ~ dat1$RSx, data = dat1)
fit2
ggsurvplot(fit2, data = dat1, pval = T,
           legend.labs = c("l1x ", "l2x"), 
           xlab = "time(months)", 
           ylab = "surviVal of risk score",
           legend.title = "Risk score")























