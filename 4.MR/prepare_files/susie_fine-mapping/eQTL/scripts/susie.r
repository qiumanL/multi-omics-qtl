library(susieR)
library(dplyr)
library(ggplot2)
packageVersion('susieR')
set.seed(1)
#attach(N3finemapping)

args<-commandArgs()
gene<-args[6]
qtl<-args[7]
#head(N3finemapping$X)
#print(sets)
#print(sets)
sumstats.tmp<-read.table(paste0(gene,'.',qtl),header=F,sep="\t",check.names=F)
colnames(sumstats.tmp)<-c('gene','snp','betahat','sebetahat')
snp_remove<-sumstats.tmp[which(sumstats.tmp$sebetahat==0),]$snp
sumstats<-sumstats.tmp[!sumstats.tmp$snp %in% snp_remove,]
dim(sumstats)
qtl.snp<-sumstats.tmp$snp
##
genotype.tmp<-read.table(paste0(gene,'.genotype'),header=T,sep=" ",row.names=1,check.names=F)
#genotype<-t(genotype)
genotype<-select(genotype.tmp,-snp_remove)
dim(genotype)
genotype.snp<-colnames(genotype)
####
all.snp<-intersect(qtl.snp,genotype.snp)
sumstats<-sumstats[sumstats$snp %in% all.snp,]
genotype<-select(genotype,all.snp)

genotype.sca<-scale(genotype)
#head(genotype.sca)
## in-sample LD
R <- cor(genotype.sca)
### check the inconsistency of LD-matrix with sumstat
#z_scores <- sumstats$betahat / sumstats$sebetahat
#Rin = R
#attr(Rin, "eigen") = eigen(Rin, symmetric = TRUE)
#lambda = estimate_s_rss(z_scores, Rin, n=185)
#lambda
#condz_in = kriging_rss(z_scores, Rin, n=185)
#ggsave('ld-matrix.diagnosis.png',condz_in$plot)
## exp
exp<-as.numeric(read.table(paste0(gene,'.exp'),sep="\t",row.names=1,check.names=F)[1,])
class(exp)
#head(exp)

##
#fitted_rss3 <- susie_rss(z_scores, R_ref, n=n, L = 10)
#head(R)
## sample number
fitted_rss1 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = 185, R = R, var_y = var(exp), L = 10,
                         estimate_residual_variance = TRUE)
#fitted_rss1
#summary(fitted_rss1)
write.table(summary(fitted_rss1)$cs,'cs.summary.txt',quote=F,row.names=F)
write.table(fitted_rss1$pip,'cs.pip.txt',quote=F,row.names=F)
