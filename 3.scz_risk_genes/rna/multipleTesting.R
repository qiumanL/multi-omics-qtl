Args <- commandArgs(TRUE)
f<-Args[1]
d<-read.csv(f,head=T)
library("qvalue")
d$qvalue<-qvalue(d$pvalue)$qvalues
d$FDR<-p.adjust(d$pvalue, method="fdr")
d$bonferroni<-p.adjust(d$pvalue, method="bonferroni")
write.table(d,file=paste(f,".multiTest",sep=""),sep="\t",row.names=F,quote=F,col.names=T)

