library(dplyr)
library(tidyr)
args<-commandArgs()
type<-args[6]
d<-list()
for (i in c(1:50)){
  d[[i]]<-read.table(paste0('shuffle',i,'/',type,'/r2.',type,'.txt'),sep="\t",head=T)
  d[[i]]$shuffle<-rep(paste0('shuffle',i),dim(d[[i]])[1])
}
data<-as.data.frame(do.call("rbind",d))
head(data)
write.table(data,paste0('r2.',type,'.all.txt'),col.name=T,row.names=F,quote=F,sep="\t")
data<- data %>% pivot_wider(id_cols=gene,names_from=shuffle,values_from=c('r2','pvalue'))
head(data)
