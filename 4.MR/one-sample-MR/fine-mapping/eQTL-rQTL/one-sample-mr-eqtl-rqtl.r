library(OneSampleMR)
argv<-commandArgs()
chr<-argv[6]
gene<-as.character(argv[7])
genotype<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/susie_fine-mapping/eQTL/chr',chr,'/',gene,'/',gene,'.genotype'),sep=" ",header=T,row.names=1)
eqtl.susie<-paste0('X',as.character(read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/susie_fine-mapping/eQTL/chr',chr,'/',gene,'/cs.snp.txt'))[,1]))
eqtl.susie
snplist<-eqtl.susie
#snplist
genotype.susie<-genotype[,which(colnames(genotype) %in% snplist)]
eqtl<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/susie_fine-mapping/eQTL/chr',chr,'/',gene,'/',gene,'.exp'),sep="\t")
rqtl<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/susie_fine-mapping/rQTL/chr',chr,'/',gene,'/',gene,'.exp'),sep="\t")
eqtl<-as.numeric(eqtl[,-1])
rqtl<-as.numeric(rqtl[,-1])
dat <- cbind(genotype.susie,eqtl, rqtl)
num.snp<-length(names(genotype.susie))
num.snp
library(ivreg)
fit1<-ivreg(as.formula(paste("rqtl ~ eqtl | ", paste(names(genotype.susie)[seq(2,num.snp, by = 1)], collapse = "+"))) , data = dat)
summary(fit1)
output_test<-summary(fit1)$diagnostics
output_fit<-summary(fit1)$coefficients
write.table(output_test,'diagnostics.txt',sep="\t",quote=F,col.names=T,row.names=T)
write.table(output_fit,'slope.txt',sep="\t",quote=F,col.names=T,row.names=T)


