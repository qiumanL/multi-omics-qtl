## permute

library(qvalue)
d<-list()
for(i in c(1:22)){
    try(d[[i]]<-read.table(paste("eqtl.permute.chr",i,sep=""),head=F,sep=" "))
    print(i)
}
dd<-do.call("rbind",d)
colnames(dd)<-c("Phe_ID", "Phe_chr", "Phe_start", "Phe_end", "Strand", "N_SNPs", "distance", "SNP_ID", "SNP_chr", "SNP_start", "SNP_END", "df", "Dummy", "para1", "para2", "P.nominal", "r_squared","slope","slope_se", "P.empirical", "P.beta")
dd<-na.omit(dd)
dd_sig<-dd[dd$qvalue<0.1,]
write.table(dd, "eqtl.permute.qvalue.txt", quote=F, row.names=F, col.names=T, sep="\t")
write.table(dd_sig, "eqtl.permute.qvalue.sig.txt", quote=F, row.names=F, col.names=T, sep="\t")



