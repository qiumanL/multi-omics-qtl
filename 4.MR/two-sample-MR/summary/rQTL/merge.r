library(dplyr)
clump<-read.table('rQTL.metaxcan.gene.clump.uvmr',sep='\t')
susie<-read.table('rQTL.metaxcan.gene.susie.uvmr',sep='\t')
data<-full_join(clump,susie,by='V1')
data[is.na(data)]='-'
colnames(data)<-c('ensemble id','ribo slope','se of ribo slope','pvalue of ribo slope','intercept','pvalue of intercept','ribo slope','se of ribo slope','pvalue of ribo slope','intercept','pvalue of intercept')
head(data)
write.table(data,'two-sample-MR-result.rQTL.txt.tmp',col.names=T,row.names=F,sep='\t',quote=F)
