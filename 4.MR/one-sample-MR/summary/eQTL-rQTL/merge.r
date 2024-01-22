library(dplyr)
clump<-read.table('metaxcan.gene.clump.eQTL-rQTL',sep='\t')
susie<-read.table('metaxcan.gene.susie.eQTL-rQTL',sep='\t')
data<-full_join(clump,susie,by='V1')
data[is.na(data)]='-'
colnames(data)<-c('ensemble id','rna-ribo slope','se of rna-ribo slope','pvalue of rna-ribo slope','intercept','pvalue of intercept','rna-ribo slope','se of rna-ribo slope','pvalue of rna-ribo slope','intercept','pvalue of intercept')
head(data)
write.table(data,'one-sample-MR-eQTL-rQTL.txt.tmp',col.names=T,row.names=F,sep='\t',quote=F)
