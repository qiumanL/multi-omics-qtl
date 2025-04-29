library(OneSampleMR)
library(ivreg)
argv<-commandArgs()
chr<-argv[6]
gene<-as.character(argv[7])


genotype<-read.table(paste0(gene,'.genotype'),sep=" ",header=T,row.names=1)

snplist<-paste0('X',as.character(read.table(paste0('1.prepare_files/plink_ld-clump/rQTL/chr',opt$chr,'/',opt$gene,'/proxy_snp'))[,1]))
snplist

genotype.final<-genotype[,which(colnames(genotype) %in% snplist)]
genotype.final<-as.data.frame(genotype.final)

rqtl<-read.table(paste0(gene,'.ribo.exp'),sep="\t",header=T,row.names=1,check.names=F)
colnames(rqtl)
pqtl<-read.table(paste0(gene,'.protein.exp'),sep="\t",header=T,row.names=1,check.names=F)
rqtl<-t(rqtl[,rownames(genotype.final)])
pqtl<-t(pqtl[,rownames(genotype.final)])

dat <- cbind(genotype.final,rqtl, pqtl)
num.snp<-length(colnames(genotype.final))
num.snp

fit1<- tsps(as.formula(paste("pqtl ~ rqtl | ", paste(colnames(genotype.final)[seq(1,num.snp, by = 1)], collapse = "+"))) , data = dat)
print(summary(fit1))
summary(fit1)$smry$coefficients["xhat",]

output_fit<-t(summary(fit1)$smry$coefficients["xhat",])
output_fit<-as.data.frame(output_fit)
output_fit$nsnps<-num.snp
write.table(output_fit,'slope.txt',sep="\t",quote=F,col.names=T,row.names=F)

### calculate instrument F-statistics
fit1<-ivreg(as.formula(paste("pqtl ~ rqtl | ", paste(colnames(genotype.final)[seq(1,num.snp, by = 1)], collapse = "+"))) , data = dat)
output_test<-summary(fit1)$diagnostics
write.table(output_test,'diagnostics.txt',sep="\t",quote=F,col.names=T,row.names=T)

