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
cs_qtl<-read.table(paste0('../',opt$gene,'.rqtl'),head=F,sep="\t")
cs<-read.table('../cs.summary.txt',head=T,sep=" ")
snplist<- c()
for(i in 1:nrow(cs)){
   cs_snp_idx<- unlist(strsplit(cs[i,]$variable,","))
   cs_snp<- cs_qtl[cs_snp_idx,]
   print(head(cs_snp[order(cs_snp$V4),]))
   cs_snp_most<- cs_snp[order(cs_snp$V4),]$V2[1]
   print(cs_snp_most)
   snplist[i]<- cs_snp_most
}
snplist


###### read our own data
gwas <-read.table(paste0('2021PGC3.sczVScontrol.sumstat2.chr',opt$chr),header=T)
gwas<-gwas[gwas$SNP %in% snplist,]
gwas_dat <- format_data(gwas,type='outcome',snp_col = "SNP",beta_col = "b",se_col = "se",effect_allele_col ="A2",other_allele_col = "A1",
 pval_col = "p",id="gwas")

rqtl<-read.table(paste0("rqtl.allP.forgwas.chr",opt$chr,".",opt$gene),sep="\t")
rqtl<-rqtl[,c(1,1,2,3,4)]
colnames(rqtl)<-c("SNP","SNP_split","p","b","se")
rqtl<-rqtl[rqtl$SNP %in% snplist,]
rqtl_tmp<-rqtl
rqtl_tmp<-cbind(rqtl_tmp,rep("rqtl",dim(rqtl_tmp)[1]));
colnames(rqtl_tmp)<-c("SNP","SNP_split","p","b","se","rqtl")
rqtl_tmp<- rqtl_tmp %>% separate(SNP_split, c("chr","pos","A1", "A2"), "_")
rqtl_tmp<- rqtl_tmp[,c(1,4,5,6,7,8,9)]
colnames(rqtl)<-c("SNP","A1","A2","freq","p","b","se","rqtl")
rqtl_dat <- format_data(rqtl_tmp,type='exposure',snp_col = "SNP",beta_col = "b",se_col = "se",
  effect_allele_col ="A2",other_allele_col = "A1",
  pval_col = "p",id="rqtl")

## harmonisation
dat <- harmonise_data(exposure_dat=rqtl_dat,outcome_dat=gwas_dat)

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
