library(TwoSampleMR)
library(dplyr);library(tidyr)
library(ggplot2);library(ieugwasr)
library(getopt)
spec <- matrix(
  c("chr", "c", 2, "integer","This is chromosome",
  "gene",  "g", 2, "integer", "This is ensemble gene id!",
  "protein",  "p", 2, "integer", "This is ensemble protein id!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
#### LD_clumped snplist
pqtl_snplist<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/plink_ld-clump/pQTL/chr',opt$chr,'/',opt$gene,'/',opt$protein,'/proxy_snp.snpid'))[,2]
#snplist
###### read our own data
gwas <-read.table('/mnt/NAS-PD/shareData/PGC/2021PGC3.sczVScontrol.sumstat2',header=T)
gwas<-gwas[gwas$SNP %in% pqtl_snplist,]
gwas_dat <- format_data(gwas,type='outcome',snp_col = "SNP",beta_col = "b",se_col = "se",effect_allele_col ="A2",other_allele_col = "A1",
 pval_col = "p",id="gwas")
gene<-opt$gene
#pqtl<-read.table("/zs32/data-analysis/liucy_group/liangqiuman/psychENCODE/shared_",head=F,sep="\t")
#colnames(pqtl)<-c("phe_id","phe_chr","phe_from","phe_to",")
pqtl<-read.table(paste0("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/pqtl/for_mr/chr",opt$chr,"/protein/pqtl.allP.forgwas.chr",opt$chr,".",opt$gene,'.',opt$protein),sep="\t")
#pqtl<-pqtl[,c(1,8,8,12,14,15)]
pqtl<-pqtl[,c(2,2,3,4,5)]
colnames(pqtl)<-c("SNP","SNP_split","p","b","se")
pqtl<-pqtl[pqtl$SNP %in% pqtl_snplist,]
pqtl_tmp<-pqtl
#pqtl_tmp<-pqtl[pqtl$gene==gene,-1]
pqtl_tmp<-cbind(pqtl_tmp,rep("pqtl",dim(pqtl_tmp)[1]));
colnames(pqtl_tmp)<-c("SNP","SNP_split","p","b","se","pqtl")
pqtl_tmp<- pqtl_tmp %>% separate(SNP_split, c("chr","pos","A1", "A2"), "_")
pqtl_tmp<- pqtl_tmp[,c(1,4,5,6,7,8,9)]
#colnames(pqtl)<-c("SNP","A1","A2","freq","p","b","se","pqtl")
#dim(pqtl_tmp)
pqtl_dat <- format_data(pqtl_tmp,type='exposure',snp_col = "SNP",beta_col = "b",se_col = "se",
  effect_allele_col ="A2",other_allele_col = "A1",
  pval_col = "p",id="pqtl")
#head(pqtl_dat)
#exposure_dat
dat <- harmonise_data(exposure_dat=pqtl_dat,outcome_dat=gwas_dat) # 对数据进行harmonisation
#head(dat)
# twosample package
print('twosample package')
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
