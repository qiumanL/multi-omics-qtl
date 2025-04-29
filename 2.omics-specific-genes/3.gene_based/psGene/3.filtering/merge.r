library(dplyr)

xsgene<-read.table('2.with_regress_other_two_omics/output/braingvex.db.extra.sig',head=T)

### without regress other two omics
rna.r2<-read.table('esGene/1.without_regress_other_two_omics/output/braingvex.db.extra',head=T,sep="\t")
ribo.r2<-read.table('rsGene/1.without_regress_other_two_omics/output/braingvex.db.extra',head=T,sep="\t")
protein.r2<-read.table('psGene/1.without_regress_other_two_omics/output/braingvex.db.extra',head=T,sep="\t")

rna.r2<-rna.r2[,c(1,4)]
ribo.r2<-ribo.r2[,c(1,4)]
protein.r2<-protein.r2[,c(1,4)]


data<- xsgene %>% left_join(rna.r2,by='gene')  %>% left_join(ribo.r2,by='gene') %>% left_join(protein.r2,by='gene')

head(data)
dim(data)

write.table(data,'braingvex.db.extra.sig.with_other_omics.corr.r2.txt',quote=F,sep="\t",row.names=F,col.names=T)
