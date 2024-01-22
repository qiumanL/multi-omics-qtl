
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

## Ribo-seq specific meta information
meta.ribo<-read.table("ribo.synapseMasterSampleTracker.txt",sep="\t",head=T,comment.char="")
rownames(meta.ribo)<-meta.ribo$BID
colnames(meta.ribo)<-c("CaseID","BID","DoubleBID_orNo","Collection","BatchNum","PreparedBy","Lysate","rRNA_removal","LibrarySubmit","sequencingStatus_lanes","Date_releaseToBionimbus","RNASeqBID","MassTissue","index_omit","elutionVolumn","seqRes","Note","seqRun","nFiles")
meta.ribo$BatchNum<-factor(meta.ribo$BatchNum, levels=sort(unique(meta.ribo$BatchNum)))
meta.ribo$sequencingStatus_lanes<-as.character(meta.ribo$sequencingStatus_lanes)

## Ribo_BID to RNASeq_BID
ribo2rna<-c()
ribo2rna[as.character(meta.ribo$BID)]<-as.character(meta.ribo$RNASeqBID)

## Lysate cluster by hclust
Lysate<-matrix(as.numeric(as.Date(meta.ribo$Lysate,format="%m/%d/%Y")), ncol=1)
rownames(Lysate)<-as.character(as.Date(meta.ribo$Lysate,format="%m/%d/%Y"))
Lysate<-unique(Lysate)
Lysate.dist<-dist(Lysate,method="euclidean")
Lysate.hclust<-hclust(Lysate.dist, method="median")
Lysate.id<-cutree(Lysate.hclust,k=17)
pdf("hclust.Lysate.pdf",width=12,height=8)
plot(Lysate.hclust, main="Cluster Dendrogram", xlab="Labrary batch", ylab="Euclidean distance")
rect.hclust(Lysate.hclust,k=17)
dev.off()

meta.ribo$rRNA_removal<-as.character(meta.ribo$rRNA_removal)
meta.ribo$rRNA_removal[meta.ribo$rRNA_removal=="12/2017"]<-"12/01/2017"
rRNA_removal<-matrix(as.numeric(as.Date(meta.ribo$rRNA_removal,format="%m/%d/%Y")), ncol=1)
rownames(rRNA_removal)<-as.character(as.Date(meta.ribo$rRNA_removal,format="%m/%d/%Y"))
rRNA_removal<-unique(rRNA_removal)
rRNA_removal.dist<-dist(rRNA_removal,method="euclidean")
rRNA_removal.hclust<-hclust(rRNA_removal.dist, method="median")
rRNA_removal.id<-cutree(rRNA_removal.hclust,k=14)
pdf("hclust.rRNA_removal.pdf",width=12,height=8)
plot(rRNA_removal.hclust, main="Cluster Dendrogram", xlab="Labrary batch", ylab="Euclidean distance")
rect.hclust(rRNA_removal.hclust,k=14)
dev.off()

meta.ribo$LibrarySubmit<-as.character(meta.ribo$LibrarySubmit)
meta.ribo$LibrarySubmit[meta.ribo$LibrarySubmit=="10/19/2017, 12/1/2017"]<-"12/01/2017"
LibrarySubmit<-matrix(as.numeric(as.Date(meta.ribo$LibrarySubmit,format="%m/%d/%Y")), ncol=1)
rownames(LibrarySubmit)<-as.character(as.Date(meta.ribo$LibrarySubmit,format="%m/%d/%Y"))
LibrarySubmit<-unique(LibrarySubmit)
LibrarySubmit.dist<-dist(LibrarySubmit,method="euclidean")
LibrarySubmit.hclust<-hclust(LibrarySubmit.dist, method="median")
LibrarySubmit.id<-cutree(LibrarySubmit.hclust,k=13)
pdf("hclust.LibrarySubmit.pdf",width=12,height=8)
plot(LibrarySubmit.hclust, main="Cluster Dendrogram", xlab="Labrary batch", ylab="Euclidean distance")
rect.hclust(LibrarySubmit.hclust,k=13)
dev.off()

## output formated meta information
meta.ribo<-cbind(meta.ribo,"LysateClustered"=paste("B",Lysate.id[as.character(as.Date(meta.ribo$Lysate,format="%m/%d/%Y"))],sep=""),"rRNA_removalClustered"=paste("B",rRNA_removal.id[as.character(as.Date(meta.ribo$rRNA_removal,format="%m/%d/%Y"))],sep=""),"LibrarySubmitClustered"=paste("B",LibrarySubmit.id[as.character(as.Date(meta.ribo$LibrarySubmit,format="%m/%d/%Y"))],sep=""))
meta.ribo$LysateClustered<-factor(meta.ribo$LysateClustered, levels=paste("B",sort(unique(Lysate.id[as.character(as.Date(meta.ribo$Lysate,format="%m/%d/%Y"))])),sep=""))
meta.ribo$rRNA_removalClustered<-factor(meta.ribo$rRNA_removalClustered, levels=paste("B",sort(unique(rRNA_removal.id[as.character(as.Date(meta.ribo$rRNA_removal,format="%m/%d/%Y"))])),sep=""))
meta.ribo$LibrarySubmitClustered<-factor(meta.ribo$LibrarySubmitClustered, levels=paste("B",sort(unique(LibrarySubmit.id[as.character(as.Date(meta.ribo$LibrarySubmit,format="%m/%d/%Y"))])),sep=""))
write.table(meta.ribo, file="ribo.synapseMasterSampleTracker2.txt", sep="\t", row.names=F, col.names=T, quote=F)


###########################################
##-------- Read Gencode          --------##
###########################################

d<-read.table("/zs32/data-analysis/liucy_group/jiangyi/database/gencode.v19/gene2coor",sep="\t",head=F)
gene2chr<-c()
gene2chr[as.character(d$V7)]<-as.character(d$V1)
gene2type<-c()
gene2type[as.character(d$V7)]<-as.character(d$V5)


#############################################
#-------- Read Ribo-seq        -------------#
#############################################

## Read gene expression level
d<-read.csv("14Mar19_featureCounts.csv",head=T,row.names=1)
ribo<-d[,2:ncol(d)]
library(stringr)
rownames(ribo)<-str_split(rownames(ribo),"\\.",simplify=TRUE)[,1]
colnames(ribo)<-sub("\\.","-",sub("^X","",colnames(ribo)))

## log2cpm transform
library(limma)
log2cpm<-voom(ribo)$E

## Classify samples by Diagnosis
coln<-colnames(log2cpm)
log2cpm.ctrl<-log2cpm[,coln[coln %in% subset(meta,meta$Diagnosis=="Control")$BID]]
log2cpm.bp<-log2cpm[,coln[coln %in% subset(meta,meta$Diagnosis=="BP")$BID]]
log2cpm.scz<-log2cpm[,coln[coln %in% subset(meta,meta$Diagnosis=="SCZ")$BID]]
write.table(log2cpm,file="log2cpm",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(log2cpm.ctrl,file="log2cpm.ctrl",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(log2cpm.bp,file="log2cpm.bp",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(log2cpm.scz,file="log2cpm.scz",sep="\t",row.names=T,col.names=NA,quote=F)

genes_to_keep = apply(log2cpm>=log2(1),1,sum) >= round(0.25 * ncol(log2cpm))  ## stringent filter: still result in a non-normal distribution
log2cpm.fgene = log2cpm[genes_to_keep,]

## pca
library(ggfortify)
library(scales)
pdf("log2cpm.fgene.pca.pdf")
pc<-prcomp(t(log2cpm.fgene))
#autoplot(pc)
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
eigsP1<-percent(eigs[1]/sum(eigs))
eigsP2<-percent(eigs[2]/sum(eigs))
plot(pc$x[,1],pc$x[,2],pch=19,col="grey",main="PCA",xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
dev.off()

library(WGCNA)
normadj <- (0.5+0.5*bicor(log2cpm.fgene, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
k <- netsummary$Connectivity; K <- k/max(k); Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
sdout <- 3.5
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

if(FALSE){
## Heatmap
library("gplots")
Lab.palette=colorRampPalette(c('white','red'))
pdf("ribo.log2cpm.fgene.heatmap.pdf")
heatmap.2(log2cpm.fgene,main="Heatmap",col=Lab.palette,density.info="none",trace="none",dendrogram="both",scale="none",margins=c(10,6))
dev.off()
}

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
##--- Convert to true BID          ------##
###########################################

## Sample ID to True sample ID
d<-read.table("/zs32/data-analysis/liucy_group/jiangyi/psychENCODE/DRAMS/sampleRelation/relatedness.highlyRelate.v4.anchor_RNASeq.switchDirection",head=T,sep="\t")
d<-d[d[,4]!="",]
d<-d[d[,1]=="ribo",c(3,4)]
ribo2rna.trueID<-ribo2rna
ribo2rna.trueID[as.character(d[,1])]<-as.character(d[,2])
table(ribo2rna==ribo2rna.trueID)
ribo2rna.falseID<-ribo2rna.trueID[ribo2rna.trueID=="Not Found"]
ribo2rna.trueID<-ribo2rna.trueID[ribo2rna.trueID!="Not Found"]
ribo2rna.falseID
head(ribo2rna.trueID)
length(ribo2rna.trueID)
#stop("test")

write.table(ribo2rna.falseID,'ribo2rna.falseID.txt',quote=F,row.names=F)
## Convert to RNASeq BID
log2cpm.fgene.fsample.qn.RNASeqBID<-log2cpm.fgene.fsample.qn
colnames(log2cpm.fgene.fsample.qn.RNASeqBID)<-ribo2rna[colnames(log2cpm.fgene.fsample.qn.RNASeqBID)]
write.table(log2cpm.fgene.fsample.qn.RNASeqBID, file="log2cpm.fgene.fsample.qn.RNASeqBID", sep="\t", row.names=T, quote=F,col.names=NA)

## Convert to true BID
log2cpm.realign<-log2cpm[,colnames(log2cpm) %in% names(ribo2rna.trueID)]
log2cpm.realign<-log2cpm.realign[,!duplicated(ribo2rna.trueID[colnames(log2cpm.realign)])]
colnames(log2cpm.realign)<-ribo2rna.trueID[colnames(log2cpm.realign)]
write.table(log2cpm.realign, file="log2cpm.realign", sep="\t", row.names=T, quote=F,col.names=NA)

## Convert to true BID
log2cpm.fgene.fsample.qn.realign<-log2cpm.fgene.fsample.qn[,colnames(log2cpm.fgene.fsample.qn) %in% names(ribo2rna.trueID)]
log2cpm.fgene.fsample.qn.realign<-log2cpm.fgene.fsample.qn.realign[,!duplicated(ribo2rna.trueID[colnames(log2cpm.fgene.fsample.qn.realign)])]
colnames(log2cpm.fgene.fsample.qn.realign)<-ribo2rna.trueID[colnames(log2cpm.fgene.fsample.qn.realign)]
write.table(log2cpm.fgene.fsample.qn.realign, file="log2cpm.fgene.fsample.qn.realign", sep="\t", row.names=T, quote=F,col.names=NA)

log2cpm.fgene.fsample.norm.qn.realign<-log2cpm.fgene.fsample.norm.qn[,colnames(log2cpm.fgene.fsample.norm.qn) %in% names(ribo2rna.trueID)]
log2cpm.fgene.fsample.norm.qn.realign<-log2cpm.fgene.fsample.norm.qn.realign[,!duplicated(ribo2rna.trueID[colnames(log2cpm.fgene.fsample.norm.qn.realign)])]
colnames(log2cpm.fgene.fsample.norm.qn.realign)<-ribo2rna.trueID[colnames(log2cpm.fgene.fsample.norm.qn.realign)]
write.table(log2cpm.fgene.fsample.norm.qn.realign, file="log2cpm.fgene.fsample.norm.qn.realign", sep="\t", row.names=T, quote=F,col.names=NA)


###########################################
##---- Estimate hidden factors        ---##
###########################################
## Questions:
## 1. Add covariates have no effect on the model. The results didn't change at all!
## 2. How to judge number of hidden factors? GTEx: 15 PEER factors were identified for tissues with fewer than 150 samples; 30 for tissues with sample sizes between 150 and 250; and 35 for tissues with more than 250 tissues. For BrainGVEX, 30 hidden factors may be good.

## Estimate hidden factors
log2cpm.fgene.fsample.qn.realign.rmXYM<-log2cpm.fgene.fsample.qn.realign[!(gene2chr[rownames(log2cpm.fgene.fsample.qn.realign)] %in% c("chrM","chrX","chrY")),]

library(peer)
expr = t(as.matrix(log2cpm.fgene.fsample.qn.realign.rmXYM))  # N rows and G columns, where N is the number of samples, and G is the number of genes. No column or row names.
dim(expr)
factorList=list()
residList=list()
meta3 = meta2[colnames(log2cpm.fgene.fsample.qn.realign.rmXYM),]

for(nFactor in c(10,20,29,30,40,50)){
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
    colnames(factors1)<-colnames(log2cpm.fgene.fsample.qn.realign.rmXYM)
    rownames(factors1)<-paste("Factor",1:nrow(factors1),sep="")
    factorList[[nFactor]]<-factors1
    #factors2<-rbind(t(meta3),factors1)
    factors2<-cbind(t(factors1),meta3)
    #write.table(t(factors2), file=paste("log2cpm.fgene.fsample.qn.realign.peerCovariates.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
    write.table(factors1, file=paste("log2cpm.fgene.fsample.qn.realign.rmXYM.peerCovariates.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
    residList[[nFactor]]<-t(residuals)
    colnames(residList[[nFactor]])<-colnames(log2cpm.fgene.fsample.qn.realign.rmXYM)
    rownames(residList[[nFactor]])<-rownames(log2cpm.fgene.fsample.qn.realign.rmXYM)
    write.table(residList[[nFactor]], file=paste("log2cpm.fgene.fsample.qn.realign.rmXYM.peerResid.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
}


factors1<-factorList[[29]]
factors1[,colnames(log2cpm.fgene.fsample.qn.realign)]
z=apply(log2cpm.fgene.fsample.qn.realign,1,function(x){residuals(lm(x~t(factors1[,colnames(log2cpm.fgene.fsample.qn.realign)])))})
log2cpm.fgene.fsample.qn.realign.lmPEER29<-t(z)+rowMeans(log2cpm.fgene.fsample.qn.realign)
write.table(log2cpm.fgene.fsample.qn.realign.lmPEER29, file="log2cpm.fgene.fsample.qn.realign.lmPEER29", sep="\t", row.names=T, quote=F,col.names=NA)

## Normalization across samples (x - mean / sd)
rMeans<-rowMeans(log2cpm.fgene.fsample.qn.realign.lmPEER29)
rSDs<-apply(log2cpm.fgene.fsample.qn.realign.lmPEER29,1,sd)
log2cpm.fgene.fsample.qn.realign.lmPEER29.norm<-apply(log2cpm.fgene.fsample.qn.realign.lmPEER29,2,function(x) (x-rMeans)/rSDs)
write.table(log2cpm.fgene.fsample.qn.realign.lmPEER29.norm, file="log2cpm.fgene.fsample.qn.realign.lmPEER29.norm", sep="\t", row.names=T, quote=F,col.names=NA)

## distribution
library(reshape)
library(ggplot2)
log2cpm.fgene.fsample.qn.realign.lmPEER29.melt<-melt(log2cpm.fgene.fsample.qn.realign.lmPEER29)
log2cpm.fgene.fsample.qn.realign.lmPEER29.norm.melt<-melt(log2cpm.fgene.fsample.qn.realign.lmPEER29.norm)
pdf("distribution.log2cpm.fgene.fsample.qn.realign.lmPEER29.pdf")
p<-ggplot(log2cpm.fgene.fsample.qn.realign.lmPEER29.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()
pdf("distribution.log2cpm.fgene.fsample.qn.realign.lmPEER29.norm.pdf")
p<-ggplot(log2cpm.fgene.fsample.qn.realign.lmPEER29.norm.melt)+geom_density(aes(x=value,group=X2),color=rgb(0.2,0.2,0.2,0.2))+theme(legend.position = "none")
print(p)
dev.off()

##################################################
##- Bayesian Information Criterion (BIC) score -##
##- Akaike   Information Criterion (AIC) score -##
##################################################
## Question: 
## 1. BIC check each covariates independently, can't consider overlapped effect? Different with PVCA?
## 2. Why choose Diagnosis in base model? I used ~1 to build base model.
## 3. Should choose covariates improved or reduced BIC for most genes? --- choose decreased BIC
## 4. BIC or AIC? BIC is too stringent, will use AIC.
## 5. How to combine all genes?

## realign sample ID for Ribo-seq sequencing batch information
nFactor<-29
phe_ribo<-c("BatchNum", "PreparedBy", "sequencingStatus_lanes", "LysateClustered", "rRNA_removalClustered", "LibrarySubmitClustered")
meta.ribo.realign<-cbind(meta.ribo[,phe_ribo],"BID.ribo"=rownames(meta.ribo))
meta.ribo.realign$BID.rna<-ribo2rna[meta.ribo.realign$BID.ribo]
meta.ribo.realign<-meta.ribo.realign[!duplicated(meta.ribo.realign$BID.rna),]
rownames(meta.ribo.realign)<-meta.ribo.realign$BID.rna
meta.ribo.realign<-meta.ribo.realign[,1:6]

d<-cbind(t(log2cpm.fgene.fsample.qn.realign),t(factorList[[nFactor]])[colnames(log2cpm.fgene.fsample.qn.realign),],meta2[colnames(log2cpm.fgene.fsample.qn.realign),],meta.ribo.realign[colnames(log2cpm.fgene.fsample.qn.realign),])
d<-d[!is.na(d$BatchNum),]

nGene<-nrow(log2cpm.fgene.fsample.qn.realign)
nSample<-ncol(log2cpm.fgene.fsample.qn.realign)

## about two days for this step
covBIC<-c()
for(i in 1:nGene){
    capture.output(resBIC<-step(lm(d[,i] ~ d$Factor1 + d$Factor2 + d$Factor3 + d$Factor4 + d$Factor5 + d$Factor6 + d$Factor7 + d$Factor8 + d$Factor9 + d$Factor10 + d$Factor11 + d$Factor12 + d$Factor13 + d$Factor14 + d$Factor15 + d$Factor16 + d$Factor17 + d$Factor18 + d$Factor19 + d$Factor20 + d$Factor21 + d$Factor22 + d$Factor23 + d$Factor24 + d$Factor25 + d$Factor26 + d$Factor27 + d$Factor28 + d$Factor29 + d$BrainBank + d$Hemisphere + d$PMI + d$BrainWeight + d$YearAutopsy + d$Sex + d$Ethnicity + d$AgeDeath + d$Diagnosis + d$PMI_square + d$BrainWeight_square + d$YearAutopsy_square + d$AgeDeath_square + d$BatchNum + d$PreparedBy + d$sequencingStatus_lanes + d$LysateClustered + d$rRNA_removalClustered + d$LibrarySubmitClustered), k=log(nSample), direction="both"), file=paste("log/log.BIC.",i,sep=""))
    covBIC<-c(covBIC, names(attr(resBIC$terms,"dataClasses"))[-1])
}

# Draw
phe3 <- c(colnames(t(factorList[[nFactor]])),colnames(meta.ribo.realign)[1:6])
pdf("ribo.covariates.bic.pdf",width=7,height=3)
props <- table(covBIC)[paste("d$", phe3, sep="")]/nGene
bp <- barplot(props,  xlab = "Covariates", ylab = "Proportion of genes decreased BIC", main="Select covariates by decreased BIC", ylim= c(0,1),col = c("blue"), las=2, cex.axis=0.5, cex.lab=0.5, cex.main=0.7, axisnames=FALSE)
#axis(1, at = bp, labels = effectsNames, xlab = "Covariates", cex.axis = 0.5, las=2)  # vertical x-axis
text(bp, par("usr")[3]-0.02, srt = 45, adj = 1,labels = phe3, xpd = TRUE,cex=0.4)  # rotate 45 x-axis
text(bp, props, labels = round(props, 3), srt = 45, pos=3, cex = 0.3) # place numbers on top of bars 
dev.off()

save.image(file="Ribo_processing.cpm.all.RData")
