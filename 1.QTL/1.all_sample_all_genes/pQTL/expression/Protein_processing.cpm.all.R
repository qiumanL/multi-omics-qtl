
options(stringsAsFactors=FALSE)
library(reshape)
library(ggplot2)

###########################################
##-------- Read meta information --------##
###########################################

## general meta information
meta<-read.table("synapse-meta-clinical-technical-data-BrainGVEX-RNAseq-TonyaCurrated-V2-GGedits3.txt",sep="\t",head=T)
rownames(meta)<-meta$BID

## phenotypes
phe<-c("YearAutopsy", "AgeDeath", "PMI", "BrainWeight", "Diagnosis", "BrainBank", "Hemisphere", "Ethnicity", "Sex")
phe2<-c("BrainBank", "Hemisphere", "Sex", "Ethnicity", "Diagnosis", "PMI", "PMI_square", "BrainWeight", "BrainWeight_square", "YearAutopsy", "YearAutopsy_square", "AgeDeath", "AgeDeath_square")
meta2 <- meta[,phe]
meta2$PMI_square <- meta2$PMI^2
meta2$BrainWeight_square <- meta2$BrainWeight^2
meta2$YearAutopsy_square <- meta2$YearAutopsy^2
meta2$AgeDeath_square <- meta2$AgeDeath^2

###########################################
##-------- Read Gencode          --------##
###########################################

d<-read.table("gene2coor",sep="\t",head=F)
gene2chr<-c()
gene2chr[as.character(d$V7)]<-as.character(d$V1)
gene2chr.hgnc<-c()
gene2chr.hgnc[as.character(d$V6)]<-as.character(d$V1)
gene2type<-c()
gene2type[as.character(d$V7)]<-as.character(d$V5)


###########################################
##-------- Read Protein expression  -----##
###########################################

d<-read.table("Final_proteomics_data_v4.txt",head=T,sep="\t",comment.char="",check.names=F,quote="")
rownames(d)<-d$UniProtID
library(limma)
log2prot<-d[,-c(1:5)]
prot2gene<-c()
prot2gene[d$UniProtID]<-d$Gene_symbol
write.table(log2prot, file="log2prot", sep="\t", row.names=T, quote=F,col.names=NA)

## distribution
library(reshape)
library(ggplot2)
prot.melt<-melt(d[,-(1:5)])
log2prot.melt<-melt(log2prot)
head(log2prot.melt)
pdf("distribution.prot.pdf")
p<-ggplot(prot.melt)+geom_density(aes(x=value,group=variable),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2prot.pdf")
p<-ggplot(log2prot.melt)+geom_density(aes(x=value,group=variable),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()

## pca
library(ggfortify)
library(scales)
pdf("log2prot.pca.pdf")
pc<-prcomp(t(log2prot))
#autoplot(pc)
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
eigsP1<-percent(eigs[1]/sum(eigs))
eigsP2<-percent(eigs[2]/sum(eigs))
plot(pc$x[,1],pc$x[,2],pch=19,col="grey",main="PCA",xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
dev.off()


###########################################
##---  filter                        ----##
###########################################

## not filtering (This study)
## keep genes with at least 1 CPM in at least 50% of the individuals (Fromer's Nature Neuroscience paper in 2016)
## keep genes with at least 10 counts in at least 50% of the individuals (Michael's Github code)
## keep genes with at least 0.1 TPM in at least 25% of the individuals (Michael's PsychENCODE Capstone1 paper)
## Keep genes with at least 0.1 RPKM in at least 10 of the individuals (GTEx)
library(DESeq2)
genes_to_keep = apply(log2prot>=0,1,sum) >= round(0.25 * ncol(log2prot))
table(genes_to_keep)  ## all proteins kept
log2prot.fgene = log2prot[genes_to_keep,]
write.table(log2prot.fgene, file="log2prot.fgene", sep="\t", row.names=T, quote=F,col.names=NA)

## Remove outlier samples (not filtering)
## Code modified from Michael Gandal's Github: https://github.com/mgandal/TSC_MIA_RNAseq/blob/master/code/step4a_Expression_Analysis.R
## Original paper: Network methods for describing sample relationships in genomic datasets: application to Huntington’s disease
library(WGCNA)
normadj <- (0.5+0.5*bicor(log2prot.fgene, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
k <- netsummary$Connectivity; K <- k/max(k); Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
sdout <- 6
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(log2prot.fgene)[outliers]); print(table(outliers))
color<-rownames(pc$x) %in% colnames(log2prot.fgene)[outliers]
color[which(color==FALSE)]<-"grey"
color[which(color==TRUE)]<-"red"
pdf("log2prot.fgene.zscore.pdf")
plot(Z.K, col = color, pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-sdout, lty=2)
dev.off()
pdf("log2prot.fgene.pca.markOutlier.pdf")
plot(pc$x[,1],pc$x[,2],pch=19,col=color,main=paste("PCA (marked Z-score < -",sdout,")",sep=""),xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
dev.off()
log2prot.fgene.fsample <- log2prot.fgene[,!outliers]  # all samples kept
write.table(log2prot.fgene.fsample, file="log2prot.fgene.fsample", sep="\t", row.names=T, quote=F,col.names=NA)


###########################################
##--- Quantile Normalization       ------##
###########################################

## QN
library(preprocessCore)
log2prot.fgene.fsample.qn<-normalize.quantiles(as.matrix(log2prot.fgene.fsample),copy=T)  # Quantile normalization across columns
rownames(log2prot.fgene.fsample.qn)<-rownames(log2prot.fgene.fsample)
colnames(log2prot.fgene.fsample.qn)<-colnames(log2prot.fgene.fsample)
write.table(log2prot.fgene.fsample.qn, file="log2prot.fgene.fsample.qn", sep="\t", row.names=T, quote=F,col.names=NA)

## pca
pdf("log2prot.fgene.fsample.qn.pca.pdf")
pc<-prcomp(t(log2prot.fgene.fsample.qn))
#autoplot(pc)
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
percVar<-eigs/sum(eigs)
eigsP1<-percent(eigs[1]/sum(eigs))
eigsP2<-percent(eigs[2]/sum(eigs))
plot(pc$x[,1],pc$x[,2],pch=19,col="grey",main="PCA (Quantile Normalized)",xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
dev.off()

## distribution
library(reshape)
library(ggplot2)
log2prot.fgene.fsample.qn.melt<-melt(log2prot.fgene.fsample.qn)
pdf("distribution.log2prot.fgene.fsample.qn.pdf")
p<-ggplot(log2prot.fgene.fsample.qn.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()



###########################################
##---- Estimate hidden factors        ---##
###########################################

## remove genes in chrX/chrY/chrM
log2prot.fgene.fsample.qn.rmXYM<-log2prot.fgene.fsample.qn[!(gene2chr.hgnc[prot2gene[rownames(log2prot.fgene.fsample.qn)]] %in% c("chrM","chrX","chrY")),]  # no gene removed.

library(peer)
expr = t(as.matrix(log2prot.fgene.fsample.qn.rmXYM))  # N rows and G columns, where N is the number of samples, and G is the number of genes. No column or row names.
dim(expr)
factorList=list()
residList=list()

for(nFactor in c(2,6,10,19,20,21,30)){
    model = PEER()  # create the model object
    PEER_setPhenoMean(model,expr)  # set the observed data
    PEER_setNk(model,nFactor) # gradient number of factors
    PEER_getNk(model)
    PEER_setAdd_mean(model, TRUE)  # include an additional factor (covariate) to account for the mean expression
    PEER_update(model)  # perform the inference
    factors = PEER_getX(model)  # inferred confounders
    weights = PEER_getW(model)  # their weights
    precision = PEER_getAlpha(model)     # precision (inverse variance) of the weights
    residuals = PEER_getResiduals(model) # the residual dataset
    #plot(precision)
    #PEER_plotModel(model)
    Variance<-c(); for(i in 2:(nFactor+1)){Variance<-c(Variance,var(factors[,i]))}; Variance<-sort(Variance, decreasing=TRUE); Variance<-100*Variance/sum(Variance)
    pdf(paste("variance-factor.factor",nFactor,".pdf",sep=""),width=5,height=5)
    plot(Variance,type="o",pch=19,col="red",xlab="Factors",main=paste("#factors: ",nFactor,sep=""))
    dev.off()
    factors1<-t(factors[,-1,drop=F])
    colnames(factors1)<-colnames(log2prot.fgene.fsample.qn.rmXYM)
    rownames(factors1)<-paste("Factor",1:nrow(factors1),sep="")
    factorList[[nFactor]]<-factors1
    write.table(factors1, file=paste("log2prot.fgene.fsample.qn.rmXYM.peerCovariates.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
    residList[[nFactor]]<-t(residuals)
    colnames(residList[[nFactor]])<-colnames(log2prot.fgene.fsample.qn.rmXYM)
    rownames(residList[[nFactor]])<-rownames(log2prot.fgene.fsample.qn.rmXYM)
    write.table(residList[[nFactor]], file=paste("log2prot.fgene.fsample.qn.rmXYM.peerResid.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
}

nFactor=19
factors1<-factorList[[nFactor]]
head(factors1)
t(factors1[,colnames(log2prot.fgene.fsample.qn),drop=F])
z=apply(log2prot.fgene.fsample.qn,1,function(x){residuals(lm(x~t(factors1[,colnames(log2prot.fgene.fsample.qn),drop=F])))})
log2prot.fgene.fsample.qn.lmPEER19<-t(z)+rowMeans(log2prot.fgene.fsample.qn)
write.table(log2prot.fgene.fsample.qn.lmPEER19, file="log2prot.fgene.fsample.qn.lmPEER19", sep="\t", row.names=T, quote=F,col.names=NA)

## Normalization across samples (x - mean / sd)
rMeans<-rowMeans(log2prot.fgene.fsample.qn.lmPEER19)
rSDs<-apply(log2prot.fgene.fsample.qn.lmPEER19,1,sd)
log2prot.fgene.fsample.qn.lmPEER19.norm<-apply(log2prot.fgene.fsample.qn.lmPEER19,2,function(x) (x-rMeans)/rSDs)
write.table(log2prot.fgene.fsample.qn.lmPEER19.norm, file="log2prot.fgene.fsample.qn.lmPEER19.norm", sep="\t", row.names=T, quote=F,col.names=NA)

library(reshape)
library(ggplot2)
log2prot.fgene.fsample.qn.lmPEER19.melt<-melt(log2prot.fgene.fsample.qn.lmPEER19)
log2prot.fgene.fsample.qn.lmPEER19.norm.melt<-melt(log2prot.fgene.fsample.qn.lmPEER19.norm)
pdf("distribution.log2prot.fgene.fsample.qn.lmPEER19.pdf")
p<-ggplot(log2prot.fgene.fsample.qn.lmPEER19.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2prot.fgene.fsample.qn.lmPEER19.norm.pdf")
p<-ggplot(log2prot.fgene.fsample.qn.lmPEER19.norm.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()

##################################################
##- Bayesian Information Criterion (BIC) score -##
##################################################
nFactor<-19
d<-cbind(t(log2prot.fgene.fsample.qn),t(factorList[[nFactor]])[colnames(log2prot.fgene.fsample.qn),],meta2[colnames(log2prot.fgene.fsample.qn),phe2])
d$PMI[is.na(d$PMI)]<-mean(d$PMI[!is.na(d$PMI)])
d$PMI_square[is.na(d$PMI_square)]<-mean(d$PMI_square[!is.na(d$PMI_square)])

nGene<-nrow(log2prot.fgene.fsample.qn)
nSample<-ncol(log2prot.fgene.fsample.qn)

## about two days for this step
covBIC<-c()
n<-1
nTest<-5000
for(i in sample(1:nGene,nTest)){
    capture.output(resBIC<-step(lm(d[,i] ~ d$Factor1 + d$Factor2 +d$Factor3 + d$Factor4 + d$Factor5 + d$Factor6 + d$Factor7 + d$Factor8 + d$Factor9 + d$Factor10 + d$Factor11 + d$Factor12 + d$Factor13 + d$Factor14 + d$Factor15 + d$Factor16 + d$Factor17 + d$Factor18 + d$Factor19 + d$BrainBank + d$Hemisphere + d$PMI + d$BrainWeight + d$YearAutopsy + d$Sex + d$Ethnicity + d$AgeDeath + d$Diagnosis + d$PMI_square + d$BrainWeight_square + d$YearAutopsy_square + d$AgeDeath_square), k=log(nSample), direction="both"), file=paste("log/log.BIC.",i,sep=""))
    covBIC<-c(covBIC, names(attr(resBIC$terms,"dataClasses"))[-1])
    print(n)
    n<-n+1
}

# Draw
phe3 <- c(colnames(t(factorList[[nFactor]])),phe2)
pdf("prot.covariates.bic.pdf",width=7,height=3)
props <- table(covBIC)[paste("d$", phe3, sep="")]/nTest
bp <- barplot(props,  xlab = "Covariates", ylab = "Proportion of genes decreased BIC", main="Select covariates by decreased BIC", ylim= c(0,1),col = c("blue"), las=2, cex.axis=0.5, cex.lab=0.5, cex.main=0.7, axisnames=FALSE)
#axis(1, at = bp, labels = effectsNames, xlab = "Covariates", cex.axis = 0.5, las=2)  # vertical x-axis
text(bp, par("usr")[3]-0.02, srt = 45, adj = 1,labels = phe3, xpd = TRUE,cex=0.4)  # rotate 45 x-axis
text(bp, props, labels = round(props, 3), srt = 45, pos=3, cex = 0.3) # place numbers on top of bars 
dev.off()


save.image("Protein_processing.cpm.all.RData")

