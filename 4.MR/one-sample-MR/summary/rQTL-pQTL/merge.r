library(dplyr)
clump<-read.table('metaxcan.gene.clump.rQTL-pQTL',sep='\t')
susie<-read.table('metaxcan.gene.susie.rQTL-pQTL',sep='\t')
data<-full_join(clump,susie,by=c('V1','V2'))
data[is.na(data)]='-'
colnames(data)<-c('ensemble id','protein id','ribo-protein slope','se of ribo-protein slope','pvalue of ribo-protein slope','intercept','pvalue of intercept','ribo-protein slope','se of ribo-protein slope','pvalue of ribo-protein slope','intercept','pvalue of intercept')
head(data)
write.table(data,'one-sample-MR-rQTL-pQTL.txt.tmp',col.names=T,row.names=F,sep='\t',quote=F)
