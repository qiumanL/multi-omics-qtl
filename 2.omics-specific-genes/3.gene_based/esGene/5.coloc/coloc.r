#### this coloc is for the remaining 629 candidate genes that potentially contain both omics-specific and shared QTL signals or contain independent omics-specific signals from multiple omics types 

library(coloc)
library(getopt)
spec <- matrix(
    c('chr','c',1,'character','chromosome',
      'gene','g',1,'character','gene',
      'cs','s',2,'integer','credible set number'),
    byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
rqtl<-read.table('rqtl.dataset',sep="\t",head=F)
gwas<-read.table('gwas.dataset',sep="\t",head=F)
colnames(rqtl)<-c("beta","varbeta","snp","position","N","MAF")
colnames(gwas)<-c("beta","varbeta","snp","position")
rqtl$type<-'quant'
gwas$type<-'cc'
rqtl<-rqtl[rqtl$varbeta!=0,]
gwas<-gwas[gwas$varbeta!=0,]

## other cs
all_cs<-read.table(paste0('4.SuSiE/',opt$chr,'/',opt$gene,'/cs.snp.txt'),header=T,sep="\t")
this_cs<-read.table(paste0('4.SuSiE/',opt$chr,'/',opt$gene,'/cs',opt$cs,'.snp.txt'),header=T,sep="\t")
other_cs<-setdiff(all_cs$snp,this_cs$snp)
other_cs

cs_topsnp<- head(this_cs[order(this_cs$pip,decreasing=T),],1)$snp
cs_topsnp

## genotype
genotype<- read.table(paste0('4.SuSiE/',opt$chr,'/',opt$gene,'/',opt$gene,'.genotype'),head=T,row.names=1,check.names=F)
genotype<-scale(t(genotype))
dim(genotype)
ld<- cor(t(genotype))^2
dim(genotype)
dim(ld)
#head(ld)

## SNPs in modest LD (squared correlation > 0.01) with the confidence set of interest
tmp<- ld[rownames(ld)==cs_topsnp,]
dim(tmp)
tmp.idx<- which(tmp>0.01)
all.related.snp<-names(tmp)[tmp.idx]
related.snp<-setdiff(all.related.snp,this_cs$snp)
head(names(tmp))
head(related.snp)
print('note')
length(all.related.snp)
length(related.snp)
length(this_cs$snp)

## excluding variants that are in linkage (squared correlation > 0.1) with any other confidence sets
other_cs.ld<- as.data.frame(ld[related.snp,other_cs])
dim(other_cs.ld)
other_cs.ld.final<- other_cs.ld[apply(other_cs.ld,1,function(x) all(x<0.1)),,drop=FALSE]
dim(other_cs.ld.final)
related.snp<-rownames(other_cs.ld.final)
length(related.snp)
head(related.snp)

## xqtl after including cs SNPs and related SNPs
dim(rqtl)
dim(gwas)
related.snp<-c(this_cs$snp,related.snp)
length(related.snp)
rqtl<- subset(rqtl,snp %in% related.snp)
gwas<- subset(gwas,snp %in% related.snp)
dim(rqtl)
dim(gwas)
rqtl<- subset(rqtl,!(snp %in% other_cs))
gwas<- subset(gwas,!(snp %in% other_cs))
dim(rqtl)
dim(gwas)
#head(rqtl)
#head(gwas)
check_dataset(rqtl,warn.minp=1e-10)


my.res <- coloc.abf(dataset1=rqtl,
                    dataset2=gwas)
coloc_snp<-my.res$results
coloc_snp<-coloc_snp[order(coloc_snp$SNP.PP.H4,decreasing=T),]
write.table(coloc_snp,'coloc_snp_PP4.txt',quote=F,col.names=T,row.names=F,sep="\t")
write.table(my.res$summary,'result',quote=F,col.names=F,row.names=T,sep="\t")
print(my.res)
