#library(qvalue)
dd<-read.table('two-sample-MR-result.pQTL.txt',header=T,sep="\t")
dd$pvalue.of.protein.slope
dd$fdr.of.protein.slope<-p.adjust(dd$pvalue.of.protein.slope,method="fdr",n=length(dd$pvalue.of.protein.slope))
dd$fdr.of.intercept<-p.adjust(dd$pvalue.of.intercept,method="fdr",n=length(dd$pvalue.of.intercept))
head(dd)
write.table(dd, "two-sample-MR-result.pQTL.fdr.txt", quote=F, row.names=F, col.names=T, sep="\t")
##
dd.sig<-dd[dd$fdr.of.protein.slope<0.05,]
head(dd.sig)
write.table(dd.sig, "two-sample-MR-result.pQTL.slope_sig.txt", quote=F, row.names=F, col.names=T, sep="\t")
##
dd.pleio_no<-dd[dd$fdr.of.intercept>0.05,]
write.table(dd.pleio_no, "two-sample-MR-result.pQTL.pleio_no.txt", quote=F, row.names=F, col.names=T, sep="\t")

##
dd.sig.pleio_no<-dd[dd$fdr.of.protein.slope<0.05 & dd$fdr.of.intercept>0.05,]
write.table(dd.sig.pleio_no, "two-sample-MR-result.pQTL.slope_sig.pleio_no.txt", quote=F, row.names=F, col.names=T, sep="\t")

dd.nosig.pleio_yes<-dd[dd$fdr.of.protein.slope>0.05 & dd$fdr.of.intercept<0.05,]
write.table(dd.nosig.pleio_yes, "two-sample-MR-result.pQTL.slope_nosig.pleio_yes.txt", quote=F, row.names=F, col.names=T, sep="\t")

dd.sig.pleio_yes<-dd[dd$fdr.of.protein.slope<0.05 & dd$fdr.of.intercept<0.05,]
write.table(dd.sig.pleio_yes, "two-sample-MR-result.pQTL.slope_sig.pleio_yes.txt", quote=F, row.names=F, col.names=T, sep="\t")

dd.nosig.pleio_no<-dd[dd$fdr.of.protein.slope>0.05 & dd$fdr.of.intercept>0.05,]
write.table(dd.nosig.pleio_no, "two-sample-MR-result.pQTL.slope_nosig.pleio_no.txt", quote=F, row.names=F, col.names=T, sep="\t")


