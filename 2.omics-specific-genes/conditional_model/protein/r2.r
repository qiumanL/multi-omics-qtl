args<-commandArgs()
type<-args[6]
data_before<-read.table(paste0('../../raw/',type,'/output/braingvex.predict_predicted_expression.txt2'),sep="\t",head=T)
data_after<-read.table('output/braingvex.predict_predicted_expression.txt2',sep="\t",head=T)
#output<-data.frame(
#   gene=list(),
#   correlation=list(),
#   pvalue=list()
#)
genelist=vector()
r2list=vector()
plist=vector()
n=0
for(i in c(2:dim(data_after)[2])){
  gene=colnames(data_after)[i]
 # print(gene)
  if (gene %in% colnames(data_before)){
    n=n+1
    tmp_before<-data_before[,gene]
    tmp_after<-data_after[,gene]
   # print(head(tmp_after))
    print(cor.test(tmp_before,tmp_after,method="pearson")$p.value)
    print(cor.test(tmp_before,tmp_after,method="pearson")$estimate^2)
    r2=cor.test(tmp_before,tmp_after,method="pearson")$estimate^2
    pvalue=cor.test(tmp_before,tmp_after,method="pearson")$p.value
    genelist[n]=gene
    r2list[n]=r2
    plist[n]=pvalue
  }
}
output<-data.frame(
   gene=genelist,
   r2=r2list,
   pvalue=plist
)
#print(genelist)
output$fdr<-p.adjust(output$pvalue,method="BH")
head(output)
write.table(output,paste0('r2.',type,'.txt'),sep="\t",col.names=T,row.names=F,quote=F) 
