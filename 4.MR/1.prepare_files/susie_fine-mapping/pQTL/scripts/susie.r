library(susieR)
library(dplyr)
library(ggplot2)
packageVersion('susieR')
set.seed(1)

args<-commandArgs()
gene<-args[6]
qtl<-args[7]

sumstats.tmp<-read.table(paste0(gene,'.',qtl),header=F,sep="\t",check.names=F)
colnames(sumstats.tmp)<-c('gene','snp','betahat','sebetahat')
snp_remove<-sumstats.tmp[which(sumstats.tmp$sebetahat==0),]$snp
snp_remove
sumstats<-sumstats.tmp[!sumstats.tmp$snp %in% snp_remove,]
dim(sumstats)
qtl.snp<-sumstats.tmp$snp

##
genotype.tmp<-read.table(paste0(gene,'.genotype'),header=T,sep=" ",row.names=1,check.names=F)
genotype<-select(genotype.tmp,-snp_remove)
dim(genotype)
genotype.snp<-colnames(genotype)

##
all.snp<-intersect(qtl.snp,genotype.snp)
sumstats<-sumstats[sumstats$snp %in% all.snp,]
genotype<-genotype[sumstats$snp,]
genotype.sca<-scale(genotype)
dim(sumstats)
dim(genotype.sca)

## in-sample LD
R <- cor(genotype.sca)

exp<-as.numeric(read.table(paste0(gene,'.exp'),sep="\t",row.names=1,check.names=F)[1,])
class(exp)
#head(exp)

n=185
fitted_rss1 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, R = R, var_y = var(exp), L = 10,
                         estimate_residual_variance = TRUE)
#fitted_rss1
data.sum<-summary(fitted_rss1)$cs
write.table(data.sum,'cs.summary.txt',quote=F,row.names=F)

data.pip<-as.data.frame(fitted_rss1$pip)
colnames(data.pip)<-'pip'
data.pip$snp<-sumstats$snp
data.pip<-data.pip[,c(2,1)]
write.table(data.pip,'cs.pip.txt',quote=F,row.names=F,sep="\t")


for(j in 1:length(data.sum$cs)){
  i<-data.sum$cs[j]
  tmp<- data.sum[data.sum$cs==i,]
  cs_snp<-as.numeric(unlist(strsplit(tmp$variable,',')))
  #print(cs_snp)
  data.cs_snp<-data.pip[cs_snp,]
  write.table(data.cs_snp,paste0('cs',i,'.snp.txt'),quote=F,sep="\t",col.names=T,row.names=F)
  if(j==1){
    all_cs_snp<-data.cs_snp
  }else{
    all_cs_snp<-rbind(all_cs_snp,data.cs_snp)
  }
}
write.table(all_cs_snp,'cs.snp.txt',quote=F,sep="\t",col.names=T,row.names=F)
