library(TwoSampleMR)
library(dplyr);library(tidyr)
library(ggplot2);library(ieugwasr)
library(getopt)
spec <- matrix(
  c("chr", "c", 2, "integer","This is chromosome",
  "gene",  "g", 2, "integer", "This is ensemble gene id!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

#### LD_clumped snplist
eqtl_snplist<-read.table(paste0('1.prepare_files/plink_ld-clump/eQTL/chr',opt$chr,'/',opt$gene,'/proxy_snp'))[,1]

###### read our own data
gwas <-read.table(paste0('2021PGC3.sczVScontrol.sumstat2.chr',opt$chr),header=T)
gwas<-gwas[gwas$SNP %in% eqtl_snplist,]
gwas_dat <- format_data(gwas,type='outcome',snp_col = "SNP",beta_col = "b",se_col = "se",effect_allele_col ="A2",other_allele_col = "A1",
 pval_col = "p",id="gwas")

eqtl<-read.table(paste0("eqtl.allP.forgwas.chr",opt$chr,".",opt$gene),sep="\t")
eqtl<-eqtl[,c(1,1,2,3,4)]
colnames(eqtl)<-c("SNP","SNP_split","p","b","se")
eqtl<-eqtl[eqtl$SNP %in% eqtl_snplist,]
eqtl_tmp<-eqtl
eqtl_tmp<-cbind(eqtl_tmp,rep("eqtl",dim(eqtl_tmp)[1]));
colnames(eqtl_tmp)<-c("SNP","SNP_split","p","b","se","eqtl")
eqtl_tmp<- eqtl_tmp %>% separate(SNP_split, c("chr","pos","A1", "A2"), "_")
eqtl_tmp<- eqtl_tmp[,c(1,4,5,6,7,8,9)]
colnames(eqtl_tmp)<-c("SNP","A1","A2","freq","p","b","se","eqtl")
eqtl_dat <- format_data(eqtl_tmp,type='exposure',snp_col = "SNP",beta_col = "b",se_col = "se",
  effect_allele_col ="A2",other_allele_col = "A1",
  pval_col = "p",id="eqtl")

## harmonisation
dat <- harmonise_data(exposure_dat=eqtl_dat,outcome_dat=gwas_dat)


res<-mr(dat)
res
write.table(res,'twosampleMR.result.txt',quote=F,sep="\t",col.names=T,row.names=F)
res_egger<-mr_egger_regression(b_exp =dat$beta.exposure, b_out = dat$beta.outcome, se_exp = dat$se.exposure,se_out = dat$se.outcome)
output<-data.frame(b="",se="",b_pvalue="",b_i="",se_i="",b_i_pvalue="",Q="",Q_df="",Q_pval="")
output$b<-res_egger$b;output$se<-res_egger$se;output$b_pvalue<-res_egger$pval;
output$b_i<-res_egger$b_i;output$se_i<-res_egger$se_i;output$b_i_pvalue<-res_egger$pval_i;
output$Q<-res_egger$Q; output$Q_df<-res_egger$Q_df; output$Q_pval<-res_egger$Q_pval
output
write.table(output,'twosampleMR.egger.result.txt',quote=F,sep="\t",col.names=T,row.names=F)
