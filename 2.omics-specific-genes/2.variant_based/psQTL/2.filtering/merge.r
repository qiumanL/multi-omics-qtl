library(dplyr)
library(ggplot2)
eqtl<-read.table('eqtl.nominal.qvalue.txt.psQTL',head=T,sep="\t")
rqtl<-read.table('rqtl.nominal.qvalue.txt.psQTL',head=T,sep="\t")
pqtl<-read.table('pqtl.nominal.qvalue.txt.psQTL',head=T,sep="\t")

eqtl<-eqtl[,c(1,8,12,14,15)]
rqtl<-rqtl[,c(1,8,12,14,15)]
pqtl<-pqtl[,c(1,8,12,14,15)]

pqtl$Phe_ID<-gsub('\\|.*','',pqtl$Phe_ID)

eqtl$omics<-'rna'
rqtl$omics<-'ribo'
pqtl$omics<-'protein'
head(pqtl)
data<-rbind(eqtl,rqtl,pqtl)
p<-ggplot(pqtl,aes(x=-log10(P.nominal)))+
   geom_histogram(bins=40)+
   geom_vline(xintercept=-log10(0.05),color='red')+
   theme_bw()
ggsave('psqtl.pqtl.nominalP.hist.png',p)

p<-ggplot(data,aes(x=-log10(P.nominal)))+
   geom_histogram(bins=40)+
   geom_vline(xintercept=-log10(0.05),color='red')+
   theme_bw()+
   facet_wrap(~omics,nrow=3)
ggsave('psqtl.nominalP.hist.png',p)

colnames(pqtl)<-c('Phe_ID','SNP_ID','P.nominal.protein','slope.protein','slope_se.protein','omics')

data<- eqtl %>% full_join(rqtl,by=c('Phe_ID','SNP_ID'),suffix=c('.rna','.ribo')) %>% full_join(pqtl,by=c('Phe_ID','SNP_ID'),suffix=c('','.protein'))
head(data)
write.table(data,'psqtl.with_other_omics.txt',quote=F,sep="\t",row.names=F,col.names=T)

