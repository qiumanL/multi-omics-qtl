
###########################################
##-------- Read meta information --------##
##-------- Cluster library batch --------##
###########################################

## meta information
meta<-read.table("synapse-meta-clinical-technical-data-BrainGVEX-RNAseq-TonyaCurrated-V2-GGedits2.txt",sep="\t",head=T)
meta$LibraryBatch<-as.character(meta$LibraryBatch)

## LibraryBatch cluster by hclust
LibraryBatch<-matrix(as.numeric(as.Date(meta$LibraryBatch)), ncol=1)
rownames(LibraryBatch)<-meta$LibraryBatch
LibraryBatch<-unique(LibraryBatch)
LibraryBatch.dist<-dist(LibraryBatch,method="euclidean")
LibraryBatch.hclust<-hclust(LibraryBatch.dist, method="median")
LibraryBatch.id<-cutree(LibraryBatch.hclust,k=18)
pdf("hclust.LibraryBatch.pdf",width=12,height=8)
plot(LibraryBatch.hclust, main="Cluster Dendrogram", xlab="Labrary batch", ylab="Euclidean distance")
rect.hclust(LibraryBatch.hclust,k=18)
dev.off()

## RNAIsolationBatch cluster by hclust
RNAIsolationBatch<-matrix(as.numeric(as.Date(meta$RNAIsolationBatch)), ncol=1)
rownames(RNAIsolationBatch)<-meta$RNAIsolationBatch
RNAIsolationBatch<-unique(RNAIsolationBatch)
RNAIsolationBatch.dist<-dist(RNAIsolationBatch,method="euclidean")
RNAIsolationBatch.hclust<-hclust(RNAIsolationBatch.dist, method="median")
RNAIsolationBatch.id<-cutree(RNAIsolationBatch.hclust,k=14)
pdf("hclust.RNAIsolationBatch.pdf",width=12,height=8)
plot(RNAIsolationBatch.hclust, main="Cluster Dendrogram", xlab="Labrary batch", ylab="Euclidean distance")
rect.hclust(RNAIsolationBatch.hclust,k=14)
dev.off()

## output formated meta information
meta$libraryBatchClustered<-rep("",nrow(meta))
for(i in 1:nrow(meta)){
    meta$libraryBatchClustered[i]<-paste("B",LibraryBatch.id[meta$LibraryBatch[i]],sep="")
}
meta$RNAIsolationBatchClustered<-rep("",nrow(meta))
for(i in 1:nrow(meta)){
    meta$RNAIsolationBatchClustered[i]<-paste("B",RNAIsolationBatch.id[meta$RNAIsolationBatch[i]],sep="")
}
rownames(meta)<-meta$BID
meta$PMI[which(is.na(meta$PMI))]<-mean(which(!is.na(meta$PMI)))
write.table(meta, file="synapse-meta-clinical-technical-data-BrainGVEX-RNAseq-TonyaCurrated-V2-GGedits3.txt", sep="\t", row.names=F, col.names=T, quote=F)

## phenotypes
phe<-c("BrainBank", "Hemisphere", "Sex", "Ethnicity", "Diagnosis", "TissueState", "RNAIsolationBatchClustered", "libraryBatchClustered", "ERCC_Added", "FlowcellBatch", "SequencingPlatform", "PMI", "BrainWeight", "YearAutopsy", "AgeDeath", "RIN")
phe2<-c("BrainBank", "Hemisphere", "Sex", "Ethnicity", "Diagnosis", "TissueState", "RNAIsolationBatchClustered", "libraryBatchClustered", "ERCC_Added", "FlowcellBatch", "SequencingPlatform", "PMI", "PMI_square", "BrainWeight", "BrainWeight_square", "YearAutopsy", "YearAutopsy_square", "AgeDeath", "AgeDeath_square", "RIN", "RIN_square")
meta2 <- meta[,phe]
meta2$PMI_square <- meta2$PMI^2
meta2$BrainWeight_square <- meta2$BrainWeight^2
meta2$YearAutopsy_square <- meta2$YearAutopsy^2
meta2$AgeDeath_square <- meta2$AgeDeath^2
meta2$RIN_square <- meta2$RIN^2
meta.rna<-meta2

## SeqPC
library(ggfortify)
library(scales)
seqcol<-c(19:90,92:95)
pc<-prcomp(meta[,seqcol])
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
print(percent(eigs/sum(eigs)))
seqpc<-pc$x
colnames(seqpc)<-paste("Seq",colnames(seqpc),sep="")
r<-matrix(nrow=ncol(seqpc),ncol=length(seqcol))
p<-matrix(nrow=ncol(seqpc),ncol=length(seqcol))
colnames(r)<-colnames(meta[,seqcol])
colnames(p)<-colnames(meta[,seqcol])
rownames(r)<-colnames(seqpc)
rownames(p)<-colnames(seqpc)
for(m in 1:ncol(seqpc)){
    for(n in 1:length(seqcol)){
        res<-cor.test(seqpc[,m],meta[,seqcol[n]],method="pearson")
        r[m,n]<-res$estimate
        p[m,n]<-res$p.value
    }
}
rownames(r)<-paste(rownames(r),"(",percent(eigs/sum(eigs)),")",sep="")
rownames(p)<-paste(rownames(p),"(",percent(eigs/sum(eigs)),")",sep="")
write.table(r^2,file="cor.seqFeature-seqPC.rSquare.xls",sep="\t",quote=F,col.names=NA)
write.table(p,file="cor.seqFeature-seqPC.pVale.xls",sep="\t",quote=F,col.names=NA)
library(gplots)
Lab.palette=colorRampPalette(c('white','red'),space="Lab")
pdf("cor.seqFeature-seqPC.rSquare.pdf",width=15,height=12)
heatmap.2(r^2,main="Correlation of sequencing QC matrix and seqPCs (R-square)",col=Lab.palette,density.info="none",trace="none",dendrogram="none",Colv=FALSE,Rowv=FALSE,scale="none",margins=c(15,6))
dev.off()
p[p==0] <- 1e-320
pdf("cor.seqFeature-seqPC.pValue.pdf",width=15,height=12)
heatmap.2(-log(p),main="Correlation of sequencing QC matrix and seqPCs (-log(Pvalue))",col=Lab.palette,density.info="none",trace="none",dendrogram="none",Colv=FALSE,Rowv=FALSE,scale="none",margins=c(15,6))
dev.off()

## Genotype PC
genopc<-read.table("merge2.WGS-PsychChip-affy.eigenvec",head=F,sep=" ",row.names=1,check.names=F)
genopc<-genopc[,-1]
colnames(genopc)<-paste("genoPC",1:20,sep="")


###########################################
##-------- Read Gencode          --------##
###########################################

d<-read.table("gene2coor",sep="\t",head=F)
gene2chr<-c()
gene2chr[as.character(d$V7)]<-as.character(d$V1)
gene2type<-c()
gene2type[as.character(d$V7)]<-as.character(d$V5)


###########################################
##--- Read RNA quantification data   ----##
##--- log2(CPM) normalization        ----##
###########################################

## data from RujiaDai <-- Michael Gandal <mgandal@gmail.com>
# RSEM_Quant.genes.counts.RData: 2001
#  -- BrainGVEX: 429
#  -- BrainSpan: 606
#  -- CMC_MSSM: 332
#  -- CMC_PENN: 95
#  -- CMC_PITT: 186
#  -- EpiGABA: 4
#  -- iPSC: 51
#  -- UCLA-ASD: 253
#  -- Yale-ASD: 45
# RSEM_Quant_Freeze2.genes.counts.RData: 973
#  -- BipSeq: 69
#  -- CMC_HBCC: 387
#  -- EpiGABA: 22
#  -- LIBD__szControl: 495

## Load Counts
load("RSEM_Quant.genes.counts.RData")
# counts1= counts
# load("RSEM_Quant_Freeze2.genes.counts.RData")
# counts2= counts
# counts = cbind(counts1, counts2)
# rm(counts1, counts2)

## Extract BrainGVEX samples
counts<-counts[,grepl("BrainGVEX_",colnames(counts))]
colnames(counts)<-gsub("BrainGVEX_","",colnames(counts))

## log2(CPM) normalization
library(limma)
log2cpm<-voom(counts)$E

write.table(log2cpm,file="log2cpm",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(counts,file="counts",sep="\t",row.names=T,col.names=NA,quote=F)


###########################################
##---  filter                        ----##
###########################################

## keep genes with at least 1 CPM in at least 25% of the individuals (This study)
## keep genes with at least 1 CPM in at least 50% of the individuals (Fromer's Nature Neuroscience paper in 2016)
## keep genes with at least 10 counts in at least 50% of the individuals (Michael's Github code)
## keep genes with at least 0.1 TPM in at least 25% of the individuals (Michael's PsychENCODE Capstone1 paper)
## Keep genes with at least 0.1 RPKM in at least 10 of the individuals (GTEx)
library(DESeq2)
genes_to_keep = apply(log2cpm>=0,1,sum) >= round(0.25 * ncol(log2cpm))
table(genes_to_keep)
log2cpm.fgene = log2cpm[genes_to_keep,]
write.table(log2cpm.fgene, file="log2cpm.fgene", sep="\t", row.names=T, quote=F,col.names=NA)

## pca
pdf("log2cpm.fgene.pca.pdf")
pc<-prcomp(t(log2cpm.fgene))
#autoplot(pc)
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
eigsP1<-percent(eigs[1]/sum(eigs))
eigsP2<-percent(eigs[2]/sum(eigs))
plot(pc$x[,1],pc$x[,2],pch=19,col="grey",main="PCA",xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
dev.off()

## Remove outlier samples (23 samples are outliers. Currently not remove samples.)
## Code modified from Michael Gandal's Github: https://github.com/mgandal/TSC_MIA_RNAseq/blob/master/code/step4a_Expression_Analysis.R
## Original paper: Network methods for describing sample relationships in genomic datasets: application to Huntington  s disease
library(WGCNA)
normadj <- (0.5+0.5*bicor(log2cpm.fgene, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
k <- netsummary$Connectivity; K <- k/max(k); Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
sdout <- 5
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(log2cpm.fgene)[outliers]); print(table(outliers))
color<-rownames(pc$x) %in% colnames(log2cpm.fgene)[outliers]
color[which(color==FALSE)]<-"grey"
color[which(color==TRUE)]<-"red"
pdf("log2cpm.fgene.zscore.pdf")
plot(Z.K, col = color, pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-sdout, lty=2)
dev.off()
pdf("log2cpm.fgene.pca.markOutlier.pdf")
plot(pc$x[,1],pc$x[,2],pch=19,col=color,main=paste("PCA (marked Z-score < -",sdout,")",sep=""),xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
dev.off()
log2cpm.fgene.fsample <- log2cpm.fgene[,!outliers]
write.table(log2cpm.fgene.fsample, file="log2cpm.fgene.fsample", sep="\t", row.names=T, quote=F,col.names=NA)


###########################################
##--- Quantile Normalization       ------##
###########################################

## QN
library(preprocessCore)
log2cpm.fgene.fsample.qn<-normalize.quantiles(log2cpm.fgene.fsample,copy=T)  # Quantile normalization across columns
rownames(log2cpm.fgene.fsample.qn)<-rownames(log2cpm.fgene.fsample)
colnames(log2cpm.fgene.fsample.qn)<-colnames(log2cpm.fgene.fsample)
write.table(log2cpm.fgene.fsample.qn, file="log2cpm.fgene.fsample.qn", sep="\t", row.names=T, quote=F,col.names=NA)


## pca
pdf("log2cpm.fgene.fsample.qn.pca.pdf")
pc<-prcomp(t(log2cpm.fgene.fsample.qn))
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
log2cpm.melt<-melt(log2cpm)
log2cpm.fgene.fsample.melt<-melt(log2cpm.fgene.fsample)
log2cpm.fgene.fsample.qn.melt<-melt(log2cpm.fgene.fsample.qn)
pdf("distribution.log2cpm.pdf")
p<-ggplot(log2cpm.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2cpm.fgene.fsample.pdf")
p<-ggplot(log2cpm.fgene.fsample.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2cpm.fgene.fsample.qn.pdf")
p<-ggplot(log2cpm.fgene.fsample.qn.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()

###########################################
##--- Realign sample ID            ------##
###########################################

old2new<-read.table("relatedness.highlyRelate.v4.anchor_RNASeq.switchDirection",head=T,sep="\t")
old2new<-old2new[old2new[,1]=="RNASeq",c(2,4)]

coln<-colnames(log2cpm.fgene.fsample.qn)
coln[coln %in% old2new[,2]]<-NA
coln[coln %in% old2new[,1]]<-as.character(old2new[,2])
log2cpm.fgene.fsample.qn.realign<-log2cpm.fgene.fsample.qn
colnames(log2cpm.fgene.fsample.qn.realign)<-coln
log2cpm.fgene.fsample.qn.realign<-log2cpm.fgene.fsample.qn.realign[,!is.na(colnames(log2cpm.fgene.fsample.qn.realign))]
write.table(log2cpm.fgene.fsample.qn.realign, file="log2cpm.fgene.fsample.qn.realign", sep="\t", row.names=T, quote=F,col.names=NA)

## PCA
pdf("log2cpm.fgene.fsample.qn.realign.pca.pdf")
pc<-prcomp(t(log2cpm.fgene.fsample.qn.realign))
#autoplot(pc)
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
percVar<-eigs/sum(eigs)
eigsP1<-percent(eigs[1]/sum(eigs))
eigsP2<-percent(eigs[2]/sum(eigs))
plot(pc$x[,1],pc$x[,2],pch=19,col="grey",main="PCA",xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
dev.off()
write.table(t(pc$x),file="log2cpm.fgene.fsample.qn.realign.pca",sep="\t",row.names=T,col.names=NA,quote=F)


###########################################
##---- Estimate hidden factors        ---##
###########################################

## remove genes in chrX/chrY/chrM
log2cpm.fgene.fsample.qn.realign.rmXYM<-log2cpm.fgene.fsample.qn.realign[!(gene2chr[rownames(log2cpm.fgene.fsample.qn.realign)] %in% c("chrM","chrX","chrY")),]
samples<-colnames(log2cpm.fgene.fsample.qn.realign.rmXYM)

## Known covariates + Sequencing QC matrix PC (top 7) + genotyping PC (top 5)
meta2<-meta[samples,phe]
seqpc2<-seqpc[samples,1:7]
genopc2<-genopc[samples,1:5]
genopc2.na<-genopc2[apply(genopc2,1,function(x) any(is.na(x))),]
genopc2[apply(genopc2,1,function(x) any(is.na(x))),]<-matrix(rep(colMeans(genopc)[1:5],each=nrow(genopc2.na)),nrow=nrow(genopc2.na))
meta.model<-model.matrix( ~ 
    meta2[,"BrainBank"] + 
    meta2[,"Sex"] + 
    meta2[,"Ethnicity"] + 
    meta2[,"Diagnosis"] + 
    meta2[,"PMI"] + 
    meta2[,"YearAutopsy"] + 
    meta2[,"AgeDeath"] + 
    meta2[,"libraryBatchClustered"] + 
    meta2[,"RNAIsolationBatchClustered"] + 
    meta2[,"RIN"] + 
    seqpc2[,1] + 
    seqpc2[,2] + 
    seqpc2[,3] + 
    seqpc2[,4] + 
    seqpc2[,5] + 
    seqpc2[,6] + 
    seqpc2[,7] + 
    genopc2[,1] + 
    genopc2[,2] + 
    genopc2[,3] + 
    genopc2[,4] + 
    genopc2[,5])
rownames(meta.model)<-samples

library(peer)
expr = t(as.matrix(log2cpm.fgene.fsample.qn.realign.rmXYM))  # N rows and G columns, where N is the number of samples, and G is the number of genes. No column or row names.
dim(expr)
factorList=list()
residList=list()

for(nFactor in 1:10*10){
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
    factors1<-t(factors)[-1,]
    colnames(factors1)<-colnames(log2cpm.fgene.fsample.qn.realign.rmXYM)
    rownames(factors1)<-paste("Factor",1:nrow(factors1),sep="")
    factorList[[nFactor]]<-factors1
    factors2<-rbind(factorList[[nFactor]][,samples],t(meta.model)[,samples])
    write.table(factors2, file=paste("log2cpm.fgene.fsample.qn.realign.rmXYM.allCovariates.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
    write.table(factors1, file=paste("log2cpm.fgene.fsample.qn.realign.rmXYM.peerCovariates.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
    residList[[nFactor]]<-t(residuals)
    colnames(residList[[nFactor]])<-colnames(log2cpm.fgene.fsample.qn.realign.rmXYM)
    rownames(residList[[nFactor]])<-rownames(log2cpm.fgene.fsample.qn.realign.rmXYM)
    write.table(residList[[nFactor]], file=paste("log2cpm.fgene.fsample.qn.realign.rmXYM.peerResid.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
}

meta3<-rbind(factorList[[30]][,samples],t(meta.model)[,samples])
write.table(t(factors2), file=paste("log2cpm.fgene.fsample.qn.realign.rmXYM.allCovariates.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)

factors1<-factorList[[30]]
z=apply(log2cpm.fgene.fsample.qn.realign,1,function(x){residuals(lm(x~t(factors1[,colnames(log2cpm.fgene.fsample.qn.realign)])))})
log2cpm.fgene.fsample.qn.realign.lmPEER30<-t(z)+rowMeans(log2cpm.fgene.fsample.qn.realign)
write.table(log2cpm.fgene.fsample.qn.realign.lmPEER30, file="log2cpm.fgene.fsample.qn.realign.lmPEER30", sep="\t", row.names=T, quote=F,col.names=NA)

## Normalization across samples (x - mean / sd)
rMeans<-rowMeans(log2cpm.fgene.fsample.qn.realign.lmPEER30)
rSDs<-apply(log2cpm.fgene.fsample.qn.realign.lmPEER30,1,sd)
log2cpm.fgene.fsample.qn.realign.lmPEER30.norm<-apply(log2cpm.fgene.fsample.qn.realign.lmPEER30,2,function(x) (x-rMeans)/rSDs)
write.table(log2cpm.fgene.fsample.qn.realign.lmPEER30.norm, file="log2cpm.fgene.fsample.qn.realign.lmPEER30.norm", sep="\t", row.names=T, quote=F,col.names=NA)

library(reshape)
library(ggplot2)
log2cpm.fgene.fsample.qn.realign.lmPEER30.melt<-melt(log2cpm.fgene.fsample.qn.realign.lmPEER30)
log2cpm.fgene.fsample.qn.realign.lmPEER30.norm.melt<-melt(log2cpm.fgene.fsample.qn.realign.lmPEER30.norm)
pdf("distribution.log2cpm.fgene.fsample.qn.realign.lmPEER30.pdf")
p<-ggplot(log2cpm.fgene.fsample.qn.realign.lmPEER30.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2cpm.fgene.fsample.qn.realign.lmPEER30.norm.pdf")
p<-ggplot(log2cpm.fgene.fsample.qn.realign.lmPEER30.norm.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()


##################################################
##- Bayesian Information Criterion (BIC) score -##
##################################################

dim(t(log2cpm.fgene.fsample.qn.realign.rmXYM))
dim(t(factors2))
dim(meta.rna[colnames(meta3),])
dim(seqpc[colnames(meta3),1:7])
d<-cbind(t(log2cpm.fgene.fsample.qn.realign.rmXYM),t(factors2),meta.rna[colnames(meta3),],seqpc[colnames(meta3),1:7])

nGene<-nrow(log2cpm.fgene.fsample.qn.realign.rmXYM)
nSample<-nrow(meta4)

## about two days for this step
covBIC<-c()
for(i in 1:nGene){
    capture.output(resBIC<-step(lm(d[,i] ~ d$BrainBank + d$Hemisphere + d$PMI + d$BrainWeight + d$YearAutopsy + d$Sex + d$Ethnicity + d$AgeDeath + d$Diagnosis + d$TissueState + d$RNAIsolationBatchClustered + d$RIN + d$ERCC_Added + d$libraryBatchClustered + d$FlowcellBatch + d$SequencingPlatform + d$PMI_square + d$BrainWeight_square + d$YearAutopsy_square + d$AgeDeath_square + d$RIN_square + d$Factor1 + d$Factor2 + d$Factor3 + d$Factor4 + d$Factor5 + d$Factor6 + d$Factor7 + d$Factor8 + d$Factor9 + d$Factor10 + d$Factor11 + d$Factor12 + d$Factor13 + d$Factor14 + d$Factor15 + d$Factor16 + d$Factor17 + d$Factor18 + d$Factor19 + d$Factor20 + d$Factor21 + d$Factor22 + d$Factor23 + d$Factor24 + d$Factor25 + d$Factor26 + d$Factor27 + d$Factor28 + d$Factor29 + d$Factor30 + d$SeqPC1 + d$SeqPC2 + d$SeqPC3 + d$SeqPC4 + d$SeqPC5 + d$SeqPC6 + d$SeqPC7), k=log(nSample), direction="both"), file=paste("log/log.BIC.",i,sep=""))
    covBIC<-c(covBIC, names(attr(resBIC$terms,"dataClasses"))[-1])
}


# Draw
phe3 <- c(paste("Factor",1:30,sep=""), phe2, colnames(seqpc)[1:7])
pdf("covariates.bic.pdf",width=7,height=3)
props <- table(covBIC)[paste("d$", phe3, sep="")]/nGene
bp <- barplot(props,  xlab = "Covariates", ylab = "Proportion of genes decreased BIC", main="Select covariates by decreased BIC", ylim= c(0,1),col = c("blue"), las=2, cex.axis=0.5, cex.lab=0.5, cex.main=0.7, axisnames=FALSE)
#axis(1, at = bp, labels = effectsNames, xlab = "Covariates", cex.axis = 0.5, las=2)  # vertical x-axis
text(bp, par("usr")[3]-0.02, srt = 45, adj = 1,labels = phe3, xpd = TRUE,cex=0.4)  # rotate 45 x-axis
text(bp, props, labels = round(props, 3), srt = 45, pos=3, cex = 0.3) # place numbers on top of bars 
dev.off()


save.image(file="RNA_processing.cpm.all.RData")
