library(dplyr)
clump<-read.table('pQTL.metaxcan.gene.clump.uvmr',sep='\t')
susie<-read.table('pQTL.metaxcan.gene.susie.uvmr',sep='\t')
data<-full_join(clump,susie,by=c('V1','V2'))
data[is.na(data)]='-'
colnames(data)<-c('ensemble id','protein id','protein slope','se of protein slope','pvalue of protein slope','intercept','pvalue of intercept','protein slope','se of protein slope','pvalue of protein slope','intercept','pvalue of intercept')
head(data)
write.table(data,'two-sample-MR-result.pQTL.txt.tmp',col.names=T,row.names=F,sep='\t',quote=F)
