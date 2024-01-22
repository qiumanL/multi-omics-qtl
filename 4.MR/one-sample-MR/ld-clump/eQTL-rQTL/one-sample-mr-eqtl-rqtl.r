library(OneSampleMR)
argv<-commandArgs()
chr<-argv[6]
gene<-as.character(argv[7])
genotype<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/susie_fine-mapping/eQTL/chr',chr,'/',gene,'/',gene,'.genotype'),sep=" ",header=T,row.names=1)
#head(genotype)
#X1_74171786_T_G<-as.numeric(genotype[,1])
#X1_74172093_T_A<-as.numeric(genotype[,2])
#X1_74172169_G_A<-as.numeric(genotype[,3])
eqtl_snplist<-paste0('X',as.character(read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/plink_ld-clump/eQTL/chr',chr,'/',gene,'/proxy_snp.snpid'))[,2]))
#rqtl.susie<-paste0('X',as.character(read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/replicate_metaxcan/susie/chr',chr,'/',gene,'/rqtl/cs.snp.txt'))[,1]))
#pqtl.susie<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/replicate_metaxcan/susie/chr1/ENSG00000116791/pqtl/Q08257.cs.snp.txt')
#snplist<-union(eqtl.susie,rqtl.susie)
snplist<-eqtl_snplist
#snplist
genotype.susie<-genotype[,which(colnames(genotype) %in% snplist)]
##head(genotype.susie)
eqtl<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/susie_fine-mapping/eQTL/chr',chr,'/',gene,'/',gene,'.exp'),sep="\t")
rqtl<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/susie_fine-mapping/rQTL/chr',chr,'/',gene,'/',gene,'.exp'),sep="\t")
eqtl<-as.numeric(eqtl[,-1])
rqtl<-as.numeric(rqtl[,-1])
dat <- cbind(genotype.susie,eqtl, rqtl)
#library(ivreg)
#fit<-ivreg(rqtl ~ eqtl | X1_74171786_T_G+X1_74172093_T_A+X1_74172169_G_A)
#summary(fit)
#head(dat)
num.snp<-length(names(genotype.susie))
num.snp
library(ivreg)
fit1<-ivreg(as.formula(paste("rqtl ~ eqtl | ", paste(names(genotype.susie)[seq(2,num.snp, by = 1)], collapse = "+"))) , data = dat)
output_test<-summary(fit1)$diagnostics
output_fit<-summary(fit1)$coefficients
write.table(output_test,'diagnostics.txt',sep="\t",quote=F,col.names=T,row.names=T)
write.table(output_fit,'slope.txt',sep="\t",quote=F,col.names=T,row.names=T)


#tspslogitfit<-tsri(as.formula(paste("rqtl ~ eqtl | ", paste(names(genotype.susie)[seq(2,num.snp, by = 1)], collapse = "+"))) , data = dat, link = "identity")
#summary(tspslogitfit)

#tsrilogitfit <- tsri(rqtl ~ eqtl | X1_74171786_T_G+X1_74172093_T_A , data = dat, link = "identity")
#summary(tsrilogitfit)

