library(tidyr)
library(dplyr)
args<-commandArgs()
type<-args[6]
data<-read.table(paste0('r2.',type,'.all.txt'),sep="\t",head=T)
data<-data[,-3]
head(data)
new_data<-pivot_wider(data, id_cols = c('gene'), names_from = 'shuffle',values_from = 'r2')
head(new_data)
write.table(new_data,paste0('r2.',type,'.tmp'),sep="\t",col.names=T,row.names=F,quote=F)
