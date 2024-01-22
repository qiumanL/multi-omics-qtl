library(dplyr)
clump<-read.table('eQTL.metaxcan.gene.clump.uvmr',sep='\t')
susie<-read.table('eQTL.metaxcan.gene.susie.uvmr',sep='\t')
data<-full_join(clump,susie,by='V1')
data[is.na(data)]='-'
colnames(data)<-c('ensemble id','RNA slope','se of RNA slope','pvalue of RNA slope','intercept','pvalue of intercept','RNA slope','se of RNA slope','pvalue of RNA slope','intercept','pvalue of intercept')
head(data)
write.table(data,'two-sample-MR-result.eQTL.txt.tmp',col.names=T,row.names=F,sep='\t',quote=F)

