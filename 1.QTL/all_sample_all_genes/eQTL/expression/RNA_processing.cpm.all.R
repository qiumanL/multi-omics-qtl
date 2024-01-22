
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
genopc<-read.table("/zs32/data-analysis/liucy_group/jiangyi/psychENCODE/DNA/merge/pca/merge2.WGS-PsychChip-affy.eigenvec",head=F,sep=" ",row.names=1,check.names=F)
genopc<-genopc[,-1]
colnames(genopc)<-paste("genoPC",1:20,sep="")


###########################################
##-------- Read Gencode          --------##
###########################################

d<-read.table("/zs32/data-analysis/liucy_group/jiangyi/database/gencode.v19/gene2coor",sep="\t",head=F)
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

## Classify samples by Diagnosis
#log2cpm.ctrl<-log2cpm[,coln[coln %in% subset(meta,meta$Diagnosis=="Control")$BID]]
#log2cpm.bp<-log2cpm[,coln[coln %in% subset(meta,meta$Diagnosis=="BP")$BID]]
#log2cpm.scz<-log2cpm[,coln[coln %in% subset(meta,meta$Diagnosis=="SCZ")$BID]]
write.table(log2cpm,file="log2cpm",sep="\t",row.names=T,col.names=NA,quote=F)
#write.table(log2cpm.ctrl,file="log2cpm.ctrl",sep="\t",row.names=T,col.names=NA,quote=F)
#write.table(log2cpm.bp,file="log2cpm.bp",sep="\t",row.names=T,col.names=NA,quote=F)
#write.table(log2cpm.scz,file="log2cpm.scz",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(counts,file="counts",sep="\t",row.names=T,col.names=NA,quote=F)
#write.table(tpm,file="tpm",sep="\t",row.names=T,col.names=NA,quote=F)


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
#normadj <- (0.5+0.5*cor(log2cpm.fgene))^2
netsummary <- fundamentalNetworkConcepts(normadj);
k <- netsummary$Connectivity; K <- k/max(k); Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
sdout <- 5
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
# use bicor: 22 outlier samples (log2cpm-based), 24 outlier samples (counts-based), 23 outlier samples (tpm-based)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(log2cpm.fgene)[outliers]); print(table(outliers))
color<-rownames(pc$x) %in% colnames(log2cpm.fgene)[outliers]
color[which(color==FALSE)]<-"grey"
color[which(color==TRUE)]<-"red"
pdf("log2cpm.fgene.zscore.pdf")
#plot(Z.K, col = as.numeric(colData(dds)$Region), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
plot(Z.K, col = color, pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-sdout, lty=2)
dev.off()
pdf("log2cpm.fgene.pca.markOutlier.pdf")
plot(pc$x[,1],pc$x[,2],pch=19,col=color,main=paste("PCA (marked Z-score < -",sdout,")",sep=""),xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
dev.off()
log2cpm.fgene.fsample <- log2cpm.fgene[,!outliers]
write.table(log2cpm.fgene.fsample, file="log2cpm.fgene.fsample", sep="\t", row.names=T, quote=F,col.names=NA)

## Normalization across samples (x - mean / sd)
rMeans<-rowMeans(log2cpm.fgene.fsample)
rSDs<-apply(log2cpm.fgene.fsample,1,sd)
log2cpm.fgene.fsample.norm<-apply(log2cpm.fgene.fsample,2,function(x) (x-rMeans)/rSDs)
write.table(log2cpm.fgene.fsample.norm, file="log2cpm.fgene.fsample.norm", sep="\t", row.names=T, quote=F,col.names=NA)


###########################################
##--- Quantile Normalization       ------##
###########################################

## QN
library(preprocessCore)
log2cpm.fgene.fsample.qn<-normalize.quantiles(log2cpm.fgene.fsample,copy=T)  # Quantile normalization across columns
rownames(log2cpm.fgene.fsample.qn)<-rownames(log2cpm.fgene.fsample)
colnames(log2cpm.fgene.fsample.qn)<-colnames(log2cpm.fgene.fsample)
write.table(log2cpm.fgene.fsample.qn, file="log2cpm.fgene.fsample.qn", sep="\t", row.names=T, quote=F,col.names=NA)

log2cpm.fgene.fsample.norm.qn<-normalize.quantiles(log2cpm.fgene.fsample.norm,copy=T)  # Quantile normalization across columns
rownames(log2cpm.fgene.fsample.norm.qn)<-rownames(log2cpm.fgene.fsample.norm)
colnames(log2cpm.fgene.fsample.norm.qn)<-colnames(log2cpm.fgene.fsample.norm)
write.table(log2cpm.fgene.fsample.norm.qn, file="log2cpm.fgene.fsample.norm.qn", sep="\t", row.names=T, quote=F,col.names=NA)

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
log2cpm.fgene.fsample.norm.melt<-melt(log2cpm.fgene.fsample.norm)
log2cpm.fgene.fsample.qn.melt<-melt(log2cpm.fgene.fsample.qn)
log2cpm.fgene.fsample.norm.qn.melt<-melt(log2cpm.fgene.fsample.norm.qn)
pdf("distribution.log2cpm.pdf")
p<-ggplot(log2cpm.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2cpm.fgene.fsample.pdf")
p<-ggplot(log2cpm.fgene.fsample.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2cpm.fgene.fsample.norm.pdf")
p<-ggplot(log2cpm.fgene.fsample.norm.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2cpm.fgene.fsample.qn.pdf")
p<-ggplot(log2cpm.fgene.fsample.qn.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2cpm.fgene.fsample.norm.qn.pdf")
p<-ggplot(log2cpm.fgene.fsample.norm.qn.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()

###########################################
##--- Realign sample ID            ------##
###########################################

old2new<-read.table("/zs32/data-analysis/liucy_group/jiangyi/psychENCODE/DRAMS/sampleRelation/relatedness.highlyRelate.v4.anchor_RNASeq.switchDirection",head=T,sep="\t")
old2new<-old2new[old2new[,1]=="RNASeq",c(2,4)]

coln<-colnames(log2cpm.fgene.fsample.qn)
coln[coln %in% old2new[,2]]<-NA
coln[coln %in% old2new[,1]]<-as.character(old2new[,2])
log2cpm.fgene.fsample.qn.realign<-log2cpm.fgene.fsample.qn
colnames(log2cpm.fgene.fsample.qn.realign)<-coln
log2cpm.fgene.fsample.qn.realign<-log2cpm.fgene.fsample.qn.realign[,!is.na(colnames(log2cpm.fgene.fsample.qn.realign))]
#log2cpm.fgene.fsample.qn.realign<-log2cpm.fgene.fsample.qn.realign[,!duplicated(coln)]
write.table(log2cpm.fgene.fsample.qn.realign, file="log2cpm.fgene.fsample.qn.realign", sep="\t", row.names=T, quote=F,col.names=NA)

coln<-colnames(log2cpm.fgene.fsample.norm.qn)
coln[coln %in% old2new[,2]]<-NA
coln[coln %in% old2new[,1]]<-as.character(old2new[,2])
log2cpm.fgene.fsample.norm.qn.realign<-log2cpm.fgene.fsample.norm.qn
colnames(log2cpm.fgene.fsample.norm.qn.realign)<-coln
log2cpm.fgene.fsample.norm.qn.realign<-log2cpm.fgene.fsample.norm.qn.realign[,!is.na(colnames(log2cpm.fgene.fsample.norm.qn.realign))]
write.table(log2cpm.fgene.fsample.norm.qn.realign, file="log2cpm.fgene.fsample.norm.qn.realign", sep="\t", row.names=T, quote=F,col.names=NA)

coln<-colnames(counts)
coln[coln %in% old2new[,2]]<-NA
coln[coln %in% old2new[,1]]<-as.character(old2new[,2])
counts.realign<-counts
colnames(counts.realign)<-coln
counts.realign<-counts.realign[,!is.na(colnames(counts.realign))]
write.table(counts.realign,file="counts.realign",sep="\t",row.names=T,col.names=NA,quote=F)

# # Venn plot
# require(VennDiagram)
# pdf("outliers.venn.pdf",width=4,height=4)
# VD = venn.diagram(height=4000, width=4000, margin=.07, x=list(log2cpm=names(outliers_log2cpm[outliers_log2cpm==TRUE]), tpm=names(outliers_tpm[outliers_tpm==TRUE]), counts=names(outliers_counts[outliers_counts==TRUE])), filename=NULL, fill=rainbow(3))
# grid.newpage()
# grid.draw(VD)
# dev.off()

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


###############################################
##--- Filter genes by coding region         -##
###############################################

## filter genes: only keep protein-coding genes
log2cpm.fcoding<-log2cpm[gene2type[rownames(log2cpm)]=="protein_coding",]

## filter samples
library(WGCNA)
normadj <- (0.5+0.5*bicor(log2cpm.fcoding, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
k <- netsummary$Connectivity; K <- k/max(k); Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
sdout <- 5
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(table(outliers))
color<-rownames(pc$x) %in% colnames(log2cpm.fcoding)[outliers]
color[which(color==FALSE)]<-"grey"
color[which(color==TRUE)]<-"red"
pdf("log2cpm.fcoding.zscore.pdf")
plot(Z.K, col = color, pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-sdout, lty=2)
dev.off()
pdf("log2cpm.fcoding.pca.markOutlier.pdf")
plot(pc$x[,1],pc$x[,2],pch=19,col=color,main=paste("PCA (marked Z-score < -",sdout,")",sep=""),xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
dev.off()
log2cpm.fcoding.fsample <- log2cpm.fcoding[,!outliers]
write.table(log2cpm.fcoding.fsample, file="log2cpm.fcoding.fsample", sep="\t", row.names=T, quote=F,col.names=NA)

## realign sample IDs
old2new<-read.table("/zs32/data-analysis/liucy_group/jiangyi/psychENCODE/DRAMS/sampleRelation/relatedness.highlyRelate.v4.anchor_RNASeq.switchDirection",head=T,sep="\t")
old2new<-old2new[old2new[,1]=="RNASeq",c(2,4)]
coln<-colnames(log2cpm.fcoding.fsample)
coln[coln %in% old2new[,2]]<-NA
coln[coln %in% old2new[,1]]<-as.character(old2new[,2])
log2cpm.fcoding.fsample.realign<-log2cpm.fcoding.fsample
colnames(log2cpm.fcoding.fsample.realign)<-coln
log2cpm.fcoding.fsample.realign<-log2cpm.fcoding.fsample.realign[,!is.na(colnames(log2cpm.fcoding.fsample.realign))]
#log2cpm.fcoding.fsample.realign<-log2cpm.fcoding.fsample.realign[,!duplicated(coln)]
write.table(log2cpm.fcoding.fsample.realign, file="log2cpm.fcoding.fsample.realign", sep="\t", row.names=T, quote=F,col.names=NA)

## Normalization across samples (x - mean / sd)
rMeans<-rowMeans(log2cpm.fcoding.fsample.realign)
rSDs<-apply(log2cpm.fcoding.fsample.realign,1,sd)
log2cpm.fcoding.fsample.realign.norm<-apply(log2cpm.fcoding.fsample.realign,2,function(x) (x-rMeans)/rSDs)
write.table(log2cpm.fcoding.fsample.realign.norm, file="log2cpm.fcoding.fsample.realign.norm", sep="\t", row.names=T, quote=F,col.names=NA)

## QN
library(preprocessCore)
log2cpm.fcoding.fsample.realign.norm.qn<-normalize.quantiles(log2cpm.fcoding.fsample.realign.norm,copy=T)  # Quantile normalization across columns
rownames(log2cpm.fcoding.fsample.realign.norm.qn)<-rownames(log2cpm.fcoding.fsample.realign.norm)
colnames(log2cpm.fcoding.fsample.realign.norm.qn)<-colnames(log2cpm.fcoding.fsample.realign.norm)
write.table(log2cpm.fcoding.fsample.realign.norm.qn, file="log2cpm.fcoding.fsample.realign.norm.qn", sep="\t", row.names=T, quote=F,col.names=NA)

## PEER
log2cpm.fcoding.fsample.realign.norm.qn.rmXYM<-log2cpm.fcoding.fsample.realign.norm.qn[!(gene2chr[rownames(log2cpm.fcoding.fsample.realign.norm.qn)] %in% c("chrM","chrX","chrY")),]
library(peer)
expr = t(as.matrix(log2cpm.fcoding.fsample.realign.norm.qn.rmXYM))  # N rows and G columns, where N is the number of samples, and G is the number of genes. No column or row names.
dim(expr)
#meta3 = meta2[colnames(log2cpm.fcoding.fsample.realign.norm.qn.rmXYM),]
nFactor=30
model = PEER()  # create the model object
PEER_setPhenoMean(model,expr)  # set the observed data
#PEER_setNk(model,30)  # Set the hidden confounders. 
PEER_setNk(model,nFactor) # gradient number of factors
PEER_getNk(model)
PEER_setAdd_mean(model, TRUE)  # include an additional factor (covariate) to account for the mean expression
## PEER_setCovariates(model, as.matrix(meta3))  # adding covariates has no effect on the model?
## PEER_setNMax_iterations(model, 10000)  # If the model is not converged after the default 1,000 iterations, and the variance of the residuals keeps decreasing, choose a higher value of iterations, e.g., 10,000.
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
colnames(factors1)<-colnames(log2cpm.fcoding.fsample.realign.norm.qn.rmXYM)
rownames(factors1)<-paste("Factor",1:nrow(factors1),sep="")
#factors2<-cbind(meta3,t(factors1[,rownames(meta3)]),seqpc[rownames(meta3),1:7])
#write.table(t(factors2), file=paste("log2cpm.fcoding.fsample.realign.norm.qn.rmXYM.allCovariates.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
write.table(factors1, file=paste("log2cpm.fcoding.fsample.realign.norm.qn.rmXYM.peerCovariates.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
z=apply(log2cpm.fcoding.fsample.realign.norm.qn,1,function(x){residuals(lm(x~t(factors1[,colnames(log2cpm.fcoding.fsample.realign.norm.qn)])))})
log2cpm.fcoding.fsample.realign.norm.qn.lmPEER30<-t(z)+rowMeans(log2cpm.fcoding.fsample.realign.norm.qn)
write.table(log2cpm.fcoding.fsample.realign.norm.qn.lmPEER30, file="log2cpm.fcoding.fsample.realign.norm.qn.lmPEER30", sep="\t", row.names=T, quote=F,col.names=NA)

## distribution
library(reshape)
library(ggplot2)
log2cpm.fcoding.fsample.melt<-melt(log2cpm.fcoding.fsample)
log2cpm.fcoding.fsample.realign.norm.melt<-melt(log2cpm.fcoding.fsample.realign.norm)
log2cpm.fcoding.fsample.realign.norm.qn.melt<-melt(log2cpm.fcoding.fsample.realign.norm.qn)
log2cpm.fcoding.fsample.realign.norm.qn.lmPEER30.melt<-melt(log2cpm.fcoding.fsample.realign.norm.qn.lmPEER30)
pdf("distribution.log2cpm.fcoding.fsample.pdf")
p<-ggplot(log2cpm.fcoding.fsample.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2cpm.fcoding.fsample.realign.norm.pdf")
p<-ggplot(log2cpm.fcoding.fsample.realign.norm.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2cpm.fcoding.fsample.realign.norm.qn.pdf")
p<-ggplot(log2cpm.fcoding.fsample.realign.norm.qn.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2cpm.fcoding.fsample.realign.norm.qn.lmPEER30.pdf")
p<-ggplot(log2cpm.fcoding.fsample.realign.norm.qn.lmPEER30.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()


###########################################
##--------------  SVA -------------------##
###########################################
if(FALSE){
    # ComBatBatch<-as.matrix(read.table("log2cpm.fgene.fsample.qn.realign",head=T,row.names=1,sep="\t"))
    meta<-read.table("synapse-meta-clinical-technical-data-BrainGVEX-RNAseq-TonyaCurrated-V2-GGedits3.txt",head=T,row.names=1,sep="\t")

    meta<-meta[samplelist,]  # make the order of columns in expression matrix and rows in meta information consistent.

    meta$Ethnicity[meta$Ethnicity=="Unknown"]<-"AA"
    meta$Ethnicity<-factor(as.character(meta$Ethnicity),levels=levels(meta$Ethnicity)[1:5])
    meta$libraryBatchClustered<-as.factor(as.character(meta$libraryBatchClustered))
    meta$Diagnosis[meta$Diagnosis=="BP_(not_BP)"]<-"BP"
    meta$Diagnosis<-factor(as.character(meta$Diagnosis),levels=levels(meta$Diagnosis)[c(1,3,4)])

    library(sva)
    Lmbeta<-log2cpm.fgene.fsample.qn.realign
    #mod = model.matrix( ~ meta$BrainWeight + meta$RIN + meta$YearAutopsy + meta$Ethnicity + meta$SequencingPlatform + meta$TissueState + meta$Sex + meta$Hemisphere)
    mod = model.matrix( ~ meta$Ethnicity + meta$libraryBatchClustered + meta$BrainBank + meta$AgeDeath + meta$Sex)
    n.sv = num.sv(Lmbeta, mod)   # n.sv=30
    svobj = sva(Lmbeta, mod, n.sv=n.sv)
    modSv = cbind(mod, svobj$sv)
    rownames(modSv)<-colnames(Lmbeta)
    covariates<-t(modSv)[-1,]
    rownames(covariates)<-c(sub("meta[$]","",rownames(covariates))[1:27],paste("PC",1:30,sep=""))
    temp<-covariates[1:27,]
    temp[temp==1]="Y"
    temp[temp==0]="N"
    covariates<-rbind(temp,covariates[28:57,,drop=F])
    write.table(covariates, file="log2cpm.fgene.fsample.qn.realign.covariates.xls", sep="\t", row.names=T, quote=F,col.names=NA)

    ## lm
    mmage=model.matrix(~-1+modSv)  # 1Ϊ ؾ࣬ ൱  y=ax+b*1    û м 1   򲻿  ǽؾࡣ
    z=apply(Lmbeta,1,function(x){residuals(lm(x~mmage))})
    Lmbeta2<-t(z)+rowMeans(Lmbeta)
    write.table(Lmbeta2, file="log2cpm.fgene.fsample.qn.realign.sva", sep="\t", row.names=T, quote=F,col.names=NA)
}


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

#for(nFactor in 1:10*10){
for(nFactor in c(30)){
    model = PEER()  # create the model object
    PEER_setPhenoMean(model,expr)  # set the observed data
    #PEER_setNk(model,30)  # Set the hidden confounders. 
    PEER_setNk(model,nFactor) # gradient number of factors
    PEER_getNk(model)
    PEER_setAdd_mean(model, TRUE)  # include an additional factor (covariate) to account for the mean expression
    ## PEER_setCovariates(model, as.matrix(meta.model))  # adding covariates has no effect on the model?
    ## PEER_setNMax_iterations(model, 10000)  # If the model is not converged after the default 1,000 iterations, and the variance of the residuals keeps decreasing, choose a higher value of iterations, e.g., 10,000.
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


## test correlation of hidden and known covariates
# continuous vs. continuous: Spearman
# continuous vs. categorical: One-way anovar
# categorical vs. categorical: Chi-sqaure test of independence
nFactor<-30
factors1<-factorList[[nFactor]]
samplelist.rna<-colnames(log2cpm.fgene.fsample.qn.realign)
phe.continuous<-c(paste("Factor",1:nFactor,sep=""), "YearAutopsy", "AgeDeath", "PMI", "BrainWeight", "RIN")
phe.categorical<-c("Diagnosis", "BrainBank", "Hemisphere", "Ethnicity", "Sex", "TissueState", "RNAIsolationBatchClustered", "libraryBatchClustered", "ERCC_Added", "FlowcellBatch", "SequencingPlatform")
meta4<-cbind(t(factors1[,samplelist.rna]),meta[samplelist.rna,c(phe.continuous[-(1:nFactor)],phe.categorical)])
len1<-length(phe.continuous)
len2<-length(phe.categorical)
lenSum<-len1+len2
pval<-matrix(nrow=lenSum,ncol=lenSum)
rownames(pval)<-c(phe.continuous,phe.categorical)
colnames(pval)<-c(phe.continuous,phe.categorical)
for(m in 1:len1){
    for(n in 1:len1){
        if(m<n){
            p<-NA
        }else{
            res<-cor.test(meta4[,m],meta4[,n],method="spearman")
            p<-res$p.value
        }
        pval[m,n]<-p
    }
}
for(m in 1:len1){
    for(n in (len1+1):lenSum){
        p<-NA
        pval[m,n]<-p
    }
}
for(m in (len1+1):lenSum){
    for(n in 1:len1){
        res<-summary(aov(meta4[,n] ~ meta4[,m]))
        p<-res[[1]][5][[1]][1]
        pval[m,n]<-p
    }
}
for(m in (len1+1):lenSum){
    for(n in (len1+1):lenSum){
        if(m<n){
            p<-NA
        }else{
            x<-as.character(meta4[,m])
            y<-as.character(meta4[,n])
            res<-chisq.test(table(data.frame("X"=x, "Y"=y)))
            p<-res$p.value
        }
        pval[m,n]<-p
    }
}
write.table(pval,file=paste("cor.covariates.pValue.factor",nFactor,".xls",sep=""),sep="\t",quote=F,col.names=NA)
logP <- -log10(pval)
diag(logP)<-NA
library(gplots)
pdf(paste("cor.covariates.pval.factor",nFactor,".pdf",sep=""),width=10,height=10)
heatmap.2(logP,main="Correlation of covariates (-log10(P))",col=colorRampPalette(c('white','red'),space="Lab")(n=21),breaks=c(seq(0,80,5),90,100,150,200,300),density.info="none",trace="none",dendrogram="none",Colv=FALSE,Rowv=FALSE,scale="none",margins=c(12,6))
dev.off()

if(FALSE){
    ## test correlation of hidden factors and covariates
    library(gplots)
    rsqList<-list()
    pvalList<-list()
    meta4<-cbind(meta3,seqpc[rownames(meta3),1:7])
    for(nFactor in 1:10*10){
        factors1<-factorList[[nFactor]]
        rsq<-matrix(nrow=nFactor,ncol=ncol(meta4))
        pval<-matrix(nrow=nFactor,ncol=ncol(meta4))
        rownames(rsq)<-rownames(factors1)
        rownames(pval)<-rownames(factors1)
        colnames(rsq)<-colnames(meta4)
        colnames(pval)<-colnames(meta4)
        for(m in 1:nrow(factors1)){
            for(n in 1:ncol(meta4)){
                colnames(meta4)[n]
                mod<-model.matrix( ~ meta4[,n])
                res<-summary(lm(factors1[m,] ~ mod))
                r2<-res$adj.r.squared
                p<-pf(res$fstatistic[1],res$fstatistic[2],res$fstatistic[3],lower.tail=FALSE)
                rsq[m,n]<-r2
                pval[m,n]<-p
            }
        }
        
        write.table(rsq,file=paste("cor.hidden-known.rSquare.factor",nFactor,".xls",sep=""),sep="\t",quote=F,col.names=NA)
        write.table(pval,file=paste("cor.hidden-known.pValue.factor",nFactor,".xls",sep=""),sep="\t",quote=F,col.names=NA)
        Lab.palette=colorRampPalette(c('white','red'),space="Lab")
        pdf(paste("cor.hidden-known.rSquare.factor",nFactor,".pdf",sep=""))
        heatmap.2(rsq,main="Correlation of PEER factors with covariates (R2)",col=Lab.palette,density.info="none",trace="none",dendrogram="none",Colv=FALSE,Rowv=FALSE,scale="none",margins=c(12,6))
        dev.off()
        pdf(paste("cor.hidden-known.pValue.factor",nFactor,".pdf",sep=""))
        heatmap.2(-log(pval),main="Correlation of PEER factors with covariates (-log(P))",col=Lab.palette,density.info="none",trace="none",dendrogram="none",Colv=FALSE,Rowv=FALSE,scale="none",margins=c(12,6))
        dev.off()
        
        rsqList[[nFactor]]<-rsq
        pvalList[[nFactor]]<-pval
    }
}


##################################################
##- variance                                   -##
##################################################
if(FALSE){
  ## remove samples collected too early
  log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900<-log2cpm.fgene.fsample.qn.realign.rmXYM[,!(colnames(log2cpm.fgene.fsample.qn.realign.rmXYM) %in% rownames(meta[meta$YearAutopsy==1900,]))]
  log2cpm.fgene.fsample.qn.realign.rm1900<-log2cpm.fgene.fsample.qn.realign[,!(colnames(log2cpm.fgene.fsample.qn.realign) %in% rownames(meta[meta$YearAutopsy==1900,]))]
  ###
  factors1<-factorList[[30]]
  z=apply(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900,1,function(x){residuals(lm(x~t(factors1[,colnames(log2cpm.fgene.fsample.qn.realign)])))})
  log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30<-t(z)+rowMeans(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900)
  write.table(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30, file="log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30", sep="\t", row.names=T, quote=F,col.names=NA)
  ###
  ## Normalization across samples (x - mean / sd)
  rMeans<-rowMeans(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30)
  rSDs<-apply(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30,1,sd)
  log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30.norm<-apply(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30,2,function(x) (x-rMeans)/rSDs)
  write.table(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30.norm, file="log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30.norm", sep="\t", row.names=T, quote=F,col.names=NA)

  var.log2cpm.fgene.fsample<-apply(log2cpm.fgene.fsample,1,function(x) var(x))
  var.log2cpm.fgene.fsample.qn.realign<-apply(log2cpm.fgene.fsample.qn.realign,1,function(x) var(x))
  var.log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30<-apply(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30,1,function(x) var(x))

  nn<-names(head(var.log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30,n=30))

  ##################################################
  ##- PCA                                        -##
  ##################################################

  ## remove samples collected too early
  log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900<-log2cpm.fgene.fsample.qn.realign.rmXYM[,!(colnames(log2cpm.fgene.fsample.qn.realign.rmXYM) %in% rownames(meta[meta$YearAutopsy==1900,]))]
  log2cpm.fgene.fsample.qn.realign.rm1900<-log2cpm.fgene.fsample.qn.realign[,!(colnames(log2cpm.fgene.fsample.qn.realign) %in% rownames(meta[meta$YearAutopsy==1900,]))]

  ## PCA, no adjust, colored by known covariates
  library(ggplot2)
  library(grid)
  phe_tmp<-c("YearAutopsy", "AgeDeath", "PMI", "BrainWeight", "Diagnosis", "BrainBank", "Hemisphere", "Ethnicity", "Sex", "TissueState", "RNAIsolationBatchClustered", "libraryBatchClustered", "ERCC_Added", "SequencingPlatform", "RIN")
  meta5<-meta[colnames(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900),]
  pc<-prcomp(t(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900))
  #autoplot(pc)
  eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
  percVar<-eigs/sum(eigs)
  eigsP1<-percent(eigs[1]/sum(eigs))
  eigsP2<-percent(eigs[2]/sum(eigs))
  d<-cbind(pc$x[,1:2],meta5[rownames(pc$x),phe_tmp])
  pdf(paste("log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.pca.pdf",sep=""),width=20,height=14)
  pushViewport(viewport(layout=grid.layout(4,4)))
  vplayout<-function(x,y){viewport(layout.pos.row=x,layout.pos.col=y)}
  n=1
  for(i in phe_tmp){
      p<-ggplot(d) + geom_point(aes(x=PC1,y=PC2,col=d[,i]),alpha=0.7) + labs(x=paste("PC1 (",eigsP1,")",sep=""),y=paste("PC2 (",eigsP2,")",sep=""),title=paste("PCA (colored by ",i,")",sep=""),color=i)
      print(p,vp=vplayout(floor((n-1)/4)+1,(n-1)%%4+1))
      n=n+1
  }
  dev.off()

  ## PCA after adjust hidden factors
  z=apply(log2cpm.fgene.fsample.qn.realign.rm1900,1,function(x){residuals(lm(x~t(factors1[,colnames(log2cpm.fgene.fsample.qn.realign.rm1900)])))})
  log2cpm.fgene.fsample.qn.realign.rm1900.lmPEER30<-t(z)+rowMeans(log2cpm.fgene.fsample.qn.realign.rm1900)
  write.table(log2cpm.fgene.fsample.qn.realign.rm1900.lmPEER30, file="log2cpm.fgene.fsample.qn.realign.rm1900.lmPEER30", sep="\t", row.names=T, quote=F,col.names=NA)

  z=apply(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900,1,function(x){residuals(lm(x~t(factors1[,colnames(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900)])))})
  log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30<-t(z)+rowMeans(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900)
  pc<-prcomp(t(log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30))
  #autoplot(pc)
  eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
  percVar<-eigs/sum(eigs)
  eigsP1<-percent(eigs[1]/sum(eigs))
  eigsP2<-percent(eigs[2]/sum(eigs))
  d<-cbind(pc$x[,1:2],meta5[rownames(pc$x),phe_tmp])
  pdf("log2cpm.fgene.fsample.qn.realign.rmXYM.rm1900.lmPEER30.pca.pdf",width=20,height=16)
  pushViewport(viewport(layout=grid.layout(4,4)))
  vplayout<-function(x,y){viewport(layout.pos.row=x,layout.pos.col=y)}
  n=1
  for(i in phe_tmp){
      p<-ggplot(d) + geom_point(aes(x=PC1,y=PC2,col=d[,i]),alpha=0.7) + labs(x=paste("PC1 (",eigsP1,")",sep=""),y=paste("PC2 (",eigsP2,")",sep=""),title="PCA (adjusted hidden factors)",color=i)
      print(p,vp=vplayout(floor((n-1)/4)+1,(n-1)%%4+1))
      n=n+1
      #plot(pc$x[,1],pc$x[,2],pch=19,col="grey",main="PCA (adjusted hidden factors)",xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
  }
  dev.off()
}

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
