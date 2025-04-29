library(TwoSampleMR)
library(MendelianRandomization)
library(dplyr);library(tidyr)
library(ggplot2);library(ieugwasr)
library(getopt)
spec <- matrix(
  c("chr", "c", 2, "integer","This is chromosome",
  "gene",  "g", 2, "integer", "This is ensemble gene id!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

#### LD_clumped snplist
pqtl_snplist<-read.table(paste0('1.prepare_files/plink_ld-clump/pQTL/chr',opt$chr,'/',opt$gene,'/proxy_snp'))[,1]

###### read our own data
gwas <-read.table(paste0('2021PGC3.sczVScontrol.sumstat2.chr',opt$chr),header=T)
gwas<-gwas[gwas$SNP %in% pqtl_snplist,]
gwas_dat <- format_data(gwas,type='outcome',snp_col = "SNP",beta_col = "b",se_col = "se",effect_allele_col ="A2",other_allele_col = "A1",
 pval_col = "p",id="gwas")


pqtl<-read.table(paste0("pqtl.allP.forgwas.chr",opt$chr,".",opt$gene),sep="\t")
colnames(pqtl)<-c("SNP","SNP_split","p","b","se")
pqtl<-pqtl[pqtl$SNP %in% pqtl_snplist,]
pqtl_tmp<-pqtl
pqtl_tmp<-cbind(pqtl_tmp,rep("pqtl",dim(pqtl_tmp)[1]));
colnames(pqtl_tmp)<-c("SNP","SNP_split","p","b","se","pqtl")
pqtl_tmp<- pqtl_tmp %>% separate(SNP_split, c("chr","pos","A1", "A2"), "_")
pqtl_tmp<- pqtl_tmp[,c(1,4,5,6,7,8,9)]
colnames(pqtl_tmp)<-c("SNP","A1","A2","freq","p","b","se","pqtl")
pqtl_dat <- format_data(pqtl_tmp,type='exposure',snp_col = "SNP",beta_col = "b",se_col = "se",
  effect_allele_col ="A2",other_allele_col = "A1",
  pval_col = "p",id="pqtl")

### harmonisation
dat <- harmonise_data(exposure_dat=pqtl_dat,outcome_dat=gwas_dat) 

## accounting correlation
genotype.tmp<-read.table(paste0(opt$gene,'.genotype'),head=T,sep=" ",row.names=1,check.names=F)
genotype<- genotype.tmp[,toupper(dat$SNP)]
genotype.sca<-scale(genotype)
snp_ld<- as.matrix(cor(genotype.sca))
dim(snp_ld)


### correlated snps 
data<-mr_input(bx=as.numeric(dat$beta.exposure),bxse=as.numeric(dat$se.exposure),by=as.numeric(dat$beta.outcome),byse=as.numeric(dat$se.outcome),corr=snp_ld)
res<-mr_ivw(data, correl = TRUE)
output<-data.frame(gene=opt$gene,slope=res$Estimate,se=res$StdError,pvalue=res$Pvalue,nsnp=res$SNPs,Q=res$Heter.Stat[1],pvalue.Q=res$Heter.Stat[2])
write.table(output,'twosampleMR.result.corr.txt',quote=F,sep="\t",col.names=T,row.names=F)
# egger
res<-mr_egger(data,correl = TRUE)
output<-data.frame(gene=opt$gene,slope=res$Estimate,se=res$StdError.Est,pvalue=res$Pvalue.Est,intercept=res$Intercept,se.intercept=res$StdError.Int,pvalue.intercept=res$Pvalue.Int,nsnp=res$SNPs,Q=res$Heter.Stat[1],pvalue.Q=res$Heter.Stat[2],I2=res$I.sq)
write.table(output,'twosampleMR_egger.result.corr.txt',quote=F,sep="\t",col.names=T,row.names=F)


