
## Usage: Rscript qvalue.R infile step output_prefix

Args <- commandArgs()

infile<-Args[6]  # infile<-"pvalList.eqtl-sqtl.0.01"
step0<-as.numeric(Args[7])  # step0 <- 0.01

library("qvalue")
d<-as.numeric(read.table(infile,head=F)$V1)
d<-d[!is.na(d)]
if(max(d)==1){
    maxd=0.9
}else{
    maxd=max(d)
}
print(maxd)
pi0<-pi0est(d)
pi1<-1-pi0$pi0
print(pi1)

