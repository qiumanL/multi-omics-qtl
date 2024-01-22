library(OneSampleMR)
argv<-commandArgs()
chr<-argv[6]
gene<-as.character(argv[7])
protein<-as.character(argv[8])
genotype<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/susie_fine-mapping/rQTL/chr',chr,'/',gene,'/',gene,'.genotype'),sep=" ",header=T,row.names=1)
rqtl_snplist<-paste0('X',as.character(read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/susie_fine-mapping/rQTL/chr',chr,'/',gene,'/proxy_snp.snpid'))[,2]))
snplist<-rqtl_snplist
genotype.susie<-genotype[,which(colnames(genotype) %in% snplist)]
rqtl<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/susie_fine-mapping/rQTL/chr',chr,'/',gene,'/',gene,'.exp'),sep="\t")
pqtl<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/susie_fine-mapping/rQTL/chr',chr,'/',gene,'/',gene,'.',protein,'.exp'),sep="\t")
rqtl<-as.numeric(rqtl[,-1])
pqtl<-as.numeric(pqtl[,-1])
dat <- cbind(genotype.susie,rqtl, pqtl)
num.snp<-length(names(genotype.susie))
num.snp
library(ivreg)
fit1<-ivreg(as.formula(paste("pqtl ~ rqtl | ", paste(names(genotype.susie)[seq(2,num.snp, by = 1)], collapse = "+"))) , data = dat)
summary(fit1)
output_test<-summary(fit1)$diagnostics
output_fit<-summary(fit1)$coefficients
write.table(output_test,paste0(protein,'.diagnostics.txt'),sep="\t",quote=F,col.names=T,row.names=T)
write.table(output_fit,paste0(protein,'.slope.txt'),sep="\t",quote=F,col.names=T,row.names=T)

