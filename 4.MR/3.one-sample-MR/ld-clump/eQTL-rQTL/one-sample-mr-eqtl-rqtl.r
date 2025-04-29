library(OneSampleMR)
library(ivreg)

argv<-commandArgs()
chr<-argv[6]
gene<-as.character(argv[7])


genotype<-read.table(paste0(gene,'.genotype'),sep=" ",header=T,row.names=1)

snplist<-paste0('X',as.character(read.table(paste0('1.prepare_files/plink_ld-clump/eQTL/chr',opt$chr,'/',opt$gene,'/proxy_snp'))[,1]))
snplist

genotype.final<-genotype[,which(colnames(genotype) %in% snplist)]
genotype.final<-as.data.frame(genotype.final)

rownames(genotype.final)
eqtl<-read.table(paste0(gene,'.rna.exp'),sep="\t",header=T,row.names=1,check.names=F)
colnames(eqtl)
rqtl<-read.table(paste0(gene,'.ribo.exp'),sep="\t",header=T,row.names=1,check.names=F)
eqtl<-t(eqtl[,rownames(genotype.final)])
rqtl<-t(rqtl[,rownames(genotype.final)])

dat <- cbind(genotype.final,eqtl, rqtl)
head(dat)
num.snp<-length(colnames(genotype.final))
num.snp

fit1<- tsps(as.formula(paste("rqtl ~ eqtl | ", paste(colnames(genotype.final)[seq(1,num.snp, by = 1)], collapse = "+"))) , data = dat)
print(summary(fit1))
summary(fit1)$smry$coefficients["xhat",]

output_fit<-t(summary(fit1)$smry$coefficients["xhat",])
output_fit<-as.data.frame(output_fit)
output_fit$nsnps<-num.snp
write.table(output_fit,'slope.txt',sep="\t",quote=F,col.names=T,row.names=F)

### calculate instrument F-statistics
fit1<-ivreg(as.formula(paste("rqtl ~ eqtl | ", paste(colnames(genotype.final)[seq(1,num.snp, by = 1)], collapse = "+"))) , data = dat)
output_test<-summary(fit1)$diagnostics
write.table(output_test,'diagnostics.txt',sep="\t",quote=F,col.names=T,row.names=T)
