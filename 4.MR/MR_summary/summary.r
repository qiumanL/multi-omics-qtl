## summarizing MR results
## 2 sample MR
#library(gdata)
scz.eQTL.MR.results <- read.table("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/summary/eQTL/two-sample-MR-result.eQTL.txt",sep="\t",head=T)
scz.rQTL.MR.results <- read.table("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/summary/rQTL/two-sample-MR-result.rQTL.txt",sep="\t",head=T)
scz.pQTL.MR.results <- read.table("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/summary/pQTL/two-sample-MR-result.pQTL.txt",sep="\t",head=T)

#head(scz.pQTL.MR.results)
## SuSiE vs. LD-Clumping for instrument selection
###################
## causal effect 
png('clump_vs_susie.two-sample.rna.png')
plot(scz.eQTL.MR.results$RNA.slope,scz.eQTL.MR.results$RNA.slope.1, xlab = "LD.clump", ylab = "SuSiE", main="eQTL-sczGWAS MR causal effect, RNA risk genes")
abline(0,1)
dev.off()
cor.test(scz.eQTL.MR.results$RNA.slope,as.numeric(scz.eQTL.MR.results$RNA.slope.1))
## ribo
png('clump_vs_susie.two-sample.ribo.png')
plot(scz.rQTL.MR.results$ribo.slope,scz.rQTL.MR.results$ribo.slope.1, xlab = "LD.clump", ylab = "SuSiE", main="rQTL-sczGWAS MR causal effect, ribo risk genes")
abline(0,1)
dev.off()
cor.test(scz.rQTL.MR.results$ribo.slope,as.numeric(scz.rQTL.MR.results$ribo.slope.1))
## protein
png('clump_vs_susie.two-sample.protein.png')
plot(scz.pQTL.MR.results$protein.slope,scz.pQTL.MR.results$protein.slope.1, xlab = "LD.clump", ylab = "SuSiE", main="pQTL-sczGWAS MR causal effect, protein risk genes")
abline(0,1)
dev.off()
cor.test(scz.pQTL.MR.results$protein.slope,as.numeric(scz.pQTL.MR.results$protein.slope.1))
##  all three in one plot (supplemental figure)
png('clump_vs_susie.two-sample.all.png',width=1200,height=1200)
par(mar=c(8,8,8,10))
plot(scz.eQTL.MR.results$RNA.slope,scz.eQTL.MR.results$RNA.slope.1, xlab = "LD clump", ylab = "SuSiE", main="QTL-sczGWAS 2SMR causal effect, S-PrediXcan risk genes", xlim = c(-0.1,0.1), ylim = c(-0.2,0.2), col="green", pch=19,cex=2.5,cex.axis=2,cex.lab=2.5,cex.main=3)
points(scz.rQTL.MR.results$ribo.slope,scz.rQTL.MR.results$ribo.slope.1, col="orange", pch=19,cex=2.5)
points(scz.pQTL.MR.results$protein.slope,scz.pQTL.MR.results$protein.slope.1, col="blue", pch=19,cex=2.5)
abline(0,1)
legend(x=-0.1, y=0.2,c("mRNA", "ribosome occupancy", "protein"),col = c("green","orange","blue"), lty = 1,lwd = 2, bty = "n",cex=3)
dev.off()

##################
## eQTL
scz.eQTL.MR.results$FDR.of.RNA.slope <- p.adjust(scz.eQTL.MR.results$pvalue.of.RNA.slope, method = "BH")
scz.eQTL.MR.results$TS.MR.slope <- 0
scz.eQTL.MR.results$TS.MR.slope[which(scz.eQTL.MR.results$FDR.of.RNA.slope < 0.1)] <- 1
scz.eQTL.MR.results$TS.MR.intercept <- 0
scz.eQTL.MR.results$TS.MR.intercept[which(scz.eQTL.MR.results$pvalue.of.intercept > 0.05)] <- 1
scz.eQTL.MR.results$passed.MR <- scz.eQTL.MR.results$TS.MR.slope*scz.eQTL.MR.results$TS.MR.intercept
fail.becasue.of.Egger.index <- which(scz.eQTL.MR.results$TS.MR.slope-scz.eQTL.MR.results$TS.MR.intercept == 1)
scz.eQTL.MR.results$fail.becasue.of.Egger <- 0
scz.eQTL.MR.results$fail.becasue.of.Egger[fail.becasue.of.Egger.index] <- 1
scz.eQTL.MR.results.rep<-scz.eQTL.MR.results[scz.eQTL.MR.results$passed.MR==1,]
write.table(scz.eQTL.MR.results.rep,'scz.eQTL.MR.results.rep.txt',sep="\t",quote=F,row.names=F,col.names=T)
# rQTL
scz.rQTL.MR.results$FDR.of.ribo.slope <- p.adjust(scz.rQTL.MR.results$pvalue.of.ribo.slope, method = "BH")
scz.rQTL.MR.results$TS.MR.slope <- 0
scz.rQTL.MR.results$TS.MR.slope[which(scz.rQTL.MR.results$FDR.of.ribo.slope < 0.1)] <- 1
scz.rQTL.MR.results$TS.MR.intercept <- 0
scz.rQTL.MR.results$TS.MR.intercept[which(scz.rQTL.MR.results$pvalue.of.intercept > 0.05)] <- 1
scz.rQTL.MR.results$passed.MR <- scz.rQTL.MR.results$TS.MR.slope*scz.rQTL.MR.results$TS.MR.intercept
fail.becasue.of.Egger.index <- which(scz.rQTL.MR.results$TS.MR.slope-scz.rQTL.MR.results$TS.MR.intercept == 1)
scz.rQTL.MR.results$fail.becasue.of.Egger <- 0
scz.rQTL.MR.results$fail.becasue.of.Egger[fail.becasue.of.Egger.index] <- 1
scz.rQTL.MR.results.rep<-scz.rQTL.MR.results[scz.rQTL.MR.results$passed.MR==1,]
write.table(scz.rQTL.MR.results.rep,'scz.rQTL.MR.results.rep.txt',sep="\t",quote=F,row.names=F,col.names=T)

# pQTL
scz.pQTL.MR.results$FDR.of.protein.slope <- p.adjust(scz.pQTL.MR.results$pvalue.of.protein.slope, method = "BH")
scz.pQTL.MR.results$TS.MR.slope <- 0
scz.pQTL.MR.results$TS.MR.slope[which(scz.pQTL.MR.results$FDR.of.protein.slope < 0.1)] <- 1
scz.pQTL.MR.results$TS.MR.intercept <- 0
scz.pQTL.MR.results$TS.MR.intercept[which(scz.pQTL.MR.results$pvalue.of.intercept > 0.05)] <- 1
scz.pQTL.MR.results$passed.MR <- scz.pQTL.MR.results$TS.MR.slope*scz.pQTL.MR.results$TS.MR.intercept
fail.becasue.of.Egger.index <- which(scz.pQTL.MR.results$TS.MR.slope-scz.pQTL.MR.results$TS.MR.intercept == 1)
scz.pQTL.MR.results$fail.becasue.of.Egger <- 0
scz.pQTL.MR.results$fail.becasue.of.Egger[fail.becasue.of.Egger.index] <- 1
scz.pQTL.MR.results <- scz.pQTL.MR.results[!duplicated(scz.pQTL.MR.results$ensemble.id),]
scz.pQTL.MR.results.rep<-scz.pQTL.MR.results[scz.pQTL.MR.results$passed.MR==1,]
write.table(scz.pQTL.MR.results.rep,'scz.pQTL.MR.results.rep.txt',sep="\t",quote=F,row.names=F,col.names=T)

write.table(scz.eQTL.MR.results,'scz.eQTL.MR.results.txt',sep="\t",quote=F,row.names=F,col.names=T)
write.table(scz.rQTL.MR.results,'scz.rQTL.MR.results.txt',sep="\t",quote=F,row.names=F,col.names=T)
write.table(scz.pQTL.MR.results,'scz.pQTL.MR.results.txt',sep="\t",quote=F,row.names=F,col.names=T)
## get passed/failed test counts
number.of.tests <- length(scz.eQTL.MR.results$passed.MR)+length(scz.rQTL.MR.results$passed.MR)+length(scz.pQTL.MR.results$passed.MR)
pass.count <-sum(scz.eQTL.MR.results$passed.MR)+sum(scz.rQTL.MR.results$passed.MR)+sum(scz.pQTL.MR.results$passed.MR)
fail.count <-sum(scz.eQTL.MR.results$passed.MR == 0)+sum(scz.rQTL.MR.results$passed.MR == 0)+sum(scz.pQTL.MR.results$passed.MR == 0)
fail.becasue.of.Egger.count <- sum(scz.eQTL.MR.results$fail.becasue.of.Egger)+sum(scz.rQTL.MR.results$fail.becasue.of.Egger)+sum(scz.pQTL.MR.results$fail.becasue.of.Egger)

##
## percent passed tests
print("replicated genes %:")
pass.count
pass.count/number.of.tests
## percent failed tests
print("non replicated genes %:")
fail.count
fail.count/number.of.tests
## percent failed only because of Egger out of all fails
print("non replicated genes because of Egger %:")
fail.becasue.of.Egger.count
fail.becasue.of.Egger.count/fail.count

### merge by gene ID to get gene counts
scz.xQTL.MR.results<- merge(scz.eQTL.MR.results[,c(1,2,13:16)],scz.rQTL.MR.results[,c(1,2,13:16)], by = "ensemble.id", all = T, suffixes = c(".RNA",".ribo"))
scz.xQTL.MR.results <- merge(scz.xQTL.MR.results,scz.pQTL.MR.results[,c(1,2,14:17)], by = "ensemble.id", all = T)
#head(scz.xQTL.MR.results)

#which(is.na(scz.xQTL.MR.results$passed.MR.RNA))
#scz.xQTL.MR.results$RNA <- 1
scz.xQTL.MR.results$RNA <- scz.xQTL.MR.results$passed.MR.RNA
scz.xQTL.MR.results$RNA[which(is.na(scz.xQTL.MR.results$passed.MR.RNA))] <- 0
#scz.xQTL.MR.results$ribo <- 1
scz.xQTL.MR.results$ribo <- scz.xQTL.MR.results$passed.MR.ribo
scz.xQTL.MR.results$ribo[which(is.na(scz.xQTL.MR.results$passed.MR.ribo))] <- 0
#scz.xQTL.MR.results$protein <- 1
scz.xQTL.MR.results$protein <- scz.xQTL.MR.results$passed.MR
scz.xQTL.MR.results$protein[which(is.na(scz.xQTL.MR.results$passed.MR))] <- 0
################################
scz.xQTL.MR.results$passed.MR.RNA[which(is.na(scz.xQTL.MR.results$passed.MR.RNA))] <- 0 
scz.xQTL.MR.results$passed.MR.ribo[which(is.na(scz.xQTL.MR.results$passed.MR.ribo))] <- 0 
scz.xQTL.MR.results$passed.MR[which(is.na(scz.xQTL.MR.results$passed.MR))] <- 0 

scz.xQTL.MR.results$TS.MR.replication.count <- as.numeric(scz.xQTL.MR.results$passed.MR.RNA) + as.numeric(scz.xQTL.MR.results$passed.MR.ribo) + as.numeric(scz.xQTL.MR.results$passed.MR)
scz.xQTL.MR.results$TS.MR.replication <- scz.xQTL.MR.results$TS.MR.replication.count > 0
#write.table(scz.xQTL.MR.results,'scz.xQTL.MR.results.txt',sep="\t",quote=F,col.names=T,row.names=F)
scz.xQTL.MR.results.summary <- scz.xQTL.MR.results[,c(1,17:21)]
#scz.xQTL.MR.results.summary
write.table(scz.xQTL.MR.results.summary,'scz.xQTL.MR.results.summary.txt',sep="\t",quote=F,col.names=T,row.names=F)

## number of risk genes replicated in at least one of the three 2SMR tests
print("replicated genes in at least one of the three 2SMR tests:")
sum(scz.xQTL.MR.results.summary$TS.MR.replication)
## number of risk genes
print("number of risk genes:")
length(scz.xQTL.MR.results.summary$TS.MR.replication)
#################################
## % metaXcan risk gene replicated
print("% metaXcan risk gene replicated:")
sum(scz.xQTL.MR.results.summary$TS.MR.replication)/length(scz.xQTL.MR.results.summary$TS.MR.replication)
## % metaXcan risk gene failed to replicate
print("% metaXcan risk gene failed to replicate:")
sum(!scz.xQTL.MR.results.summary$TS.MR.replication)/length(scz.xQTL.MR.results.summary$TS.MR.replication)


################## 1 sample MR
#eQTL.rQTL.MR.results <- read.xls("one-sample-MR-eqtl-rqtl.xlsx",sheet = 1, skip = 1)
eQTL.rQTL.MR.results <- read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/summary/eQTL-rQTL/one-sample-MR-eQTL-rQTL.txt',sep="\t",head=T)
png('clump_vs_susie.one-sample.rna-ribo.png')
plot(eQTL.rQTL.MR.results$rna.ribo.slope,eQTL.rQTL.MR.results$rna.ribo.slope.1, xlab = "LD.clump", ylab = "SuSiE", main="rQTL-eQTL MR causal effect, metaXcan risk genes")
abline(0,1)
dev.off()

#rQTL.pQTL.MR.results <- read.xls("one-sample-MR-rqtl-pqtl.xlsx",sheet = 1, skip = 1)
rQTL.pQTL.MR.results <- read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/summary/rQTL-pQTL/one-sample-MR-rQTL-pQTL.txt',sep="\t",head=T)
rows.to.keep <- !duplicated(rQTL.pQTL.MR.results$gene.name)
rQTL.pQTL.MR.results<- rQTL.pQTL.MR.results[rows.to.keep,]
dim(rQTL.pQTL.MR.results)
##
png('clump_vs_susie.one-sample.ribo-protein.png')
plot(rQTL.pQTL.MR.results$ribo.protein.slope,rQTL.pQTL.MR.results$ribo.protein.slope.1, xlab = "LD.clump", ylab = "SuSiE", main="pQTL-rQTL MR causal effect, metaXcan risk genes")
abline(0,1)
dev.off()

eQTL.rQTL.MR.results$FDR.of.slope <- p.adjust(eQTL.rQTL.MR.results$pvalue.of.rna.ribo.slope, method = "BH")
causal.index <- which(eQTL.rQTL.MR.results$FDR.of.slope < 0.1 & eQTL.rQTL.MR.results$pvalue.of.intercept > 0.05)
eQTL.rQTL.MR.results$pass.MR.test <- 0
eQTL.rQTL.MR.results$pass.MR.test[causal.index] <- 1
rQTL.pQTL.MR.results$FDR.of.slope <- p.adjust(rQTL.pQTL.MR.results$pvalue.of.ribo.protein.slope, method = "BH")
causal.index <- which(rQTL.pQTL.MR.results$FDR.of.slope < 0.1 & rQTL.pQTL.MR.results$pvalue.of.intercept > 0.05)
rQTL.pQTL.MR.results$pass.MR.test <- 0
rQTL.pQTL.MR.results$pass.MR.test[causal.index] <- 1

## make a table to summarize results
eQTL.rQTL.MR.results.summary <- eQTL.rQTL.MR.results[,c(1:4,14)]
rQTL.pQTL.MR.results.summary <- rQTL.pQTL.MR.results[,c(1,2,4,5,15)]
OS.MR.results.summary <- merge(eQTL.rQTL.MR.results.summary,rQTL.pQTL.MR.results.summary, by = "ensemble.id",suffixes = c(".RNAtoRIBO",".RIBOtoPROTEIN"))
scz.xQTL.MR.results.summary <- merge(scz.xQTL.MR.results.summary,OS.MR.results.summary, by="ensemble.id")
replicated.risk.genes.MR.results.summary <- scz.xQTL.MR.results.summary[scz.xQTL.MR.results.summary$TS.MR.replication == T,]
replicated.risk.genes.MR.results.summary$pass.OSMR.test.count <- replicated.risk.genes.MR.results.summary$pass.MR.test.RNAtoRIBO + replicated.risk.genes.MR.results.summary$pass.MR.test.RIBOtoPROTEIN
write.table(replicated.risk.genes.MR.results.summary,'replicated.risk.genes.MR.results.summary.txt',sep="\t",quote=F,col.names=T,row.names=F)

### subset and print gene list
## genes pass MR at both RNA->ribo and ribo->protein 
both.pass.index <- which(replicated.risk.genes.MR.results.summary$pass.OSMR.test.count==2)
print("one-sample.both.pass:")
length(both.pass.index)
both_pathways<-replicated.risk.genes.MR.results.summary[both.pass.index,]
#both_pathways
write.table(both_pathways,'both_pathways.txt',quote=F,col.names=T,row.names=F,sep="\t")
## genes pass MR at RNA->ribo but not ribo->protein 
up.pass.index <- which(replicated.risk.genes.MR.results.summary$pass.OSMR.test.count == 1 & replicated.risk.genes.MR.results.summary$pass.MR.test.RNAtoRIBO == 1)
print("one-sample.up.pass:")
length(up.pass.index)
rna_ribo_pathway<-replicated.risk.genes.MR.results.summary[up.pass.index,]
#rna_ribo_pathway
write.table(rna_ribo_pathway,'rna_ribo_pathway.txt',quote=F,col.names=T,row.names=F,sep="\t")
## genes pass MR at ribo->protein but not RNA->ribo 
down.pass.index <- which(replicated.risk.genes.MR.results.summary$pass.OSMR.test.count == 1 & replicated.risk.genes.MR.results.summary$pass.MR.test.RIBOtoPROTEIN == 1)
print("one-sample.down.pass:")
length(down.pass.index)
ribo_protein_pathway<-replicated.risk.genes.MR.results.summary[down.pass.index,]
ribo_protein_pathway
write.table(ribo_protein_pathway,'ribo_protein.pathway.txt',quote=F,col.names=T,row.names=F,sep="\t")
## genes that fail MA at both RNA->ribo and ribo->protein 
no.pass.index <- which(replicated.risk.genes.MR.results.summary$pass.OSMR.test.count==0)
print("one-sample.no.pass:")
length(no.pass.index)
no_pathway<-replicated.risk.genes.MR.results.summary[no.pass.index,]
#no_pathway
write.table(no_pathway,'no_pathway.txt',quote=F,col.names=T,row.names=F,sep="\t")

### read in instrument diagnostic results
eQTL.instrument.F <- read.table("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/summary/eQTL-rQTL/weak-instrument.diagnosis.ld-clump.txt", header = T)
eQTL.instrument.F <- eQTL.instrument.F[,c(2,5)]
rQTL.instrument.F <- read.table("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/summary/rQTL-pQTL/weak-instrument.diagnosis.ld-clump.txt", header = T)
sum(duplicated(rQTL.instrument.F$gene))
dim(rQTL.instrument.F)
rQTL.instrument.F <- rQTL.instrument.F[,c(2,5)]
instrument.F <- merge(eQTL.instrument.F,rQTL.instrument.F, by = "gene",suffixes = c(".RNA",".ribo"))
#head(instrument.F)
replicated.risk.genes.OSMR.results.summary <- merge(replicated.risk.genes.MR.results.summary,instrument.F,by.x="ensemble.id", by.y="gene")
write.table(replicated.risk.genes.OSMR.results.summary,'replicated.risk.genes.OSMR.results.summary.txt',sep="\t",col.names=T,row.names=F,quote=F)

## create index to categorize genes according to MR results 
both.pass.index <- which(replicated.risk.genes.OSMR.results.summary$pass.OSMR.test.count==2)
up.pass.index <- which(replicated.risk.genes.OSMR.results.summary$pass.OSMR.test.count == 1 & replicated.risk.genes.OSMR.results.summary$pass.MR.test.RNAtoRIBO == 1)
no.pass.index <- which(replicated.risk.genes.OSMR.results.summary$pass.OSMR.test.count==0)
### comparing ribo -> protein MR stats between both pass gene and one pass genes (note that all one pass genes pass upstream but fail downstream MR)
png('effect_size.downstream.png')
boxplot(replicated.risk.genes.OSMR.results.summary$ribo.protein.slope[both.pass.index], replicated.risk.genes.OSMR.results.summary$ribo.protein.slope[up.pass.index], names = c("pass", "fail"), main="Downstream MR: causal effect slope", outline = F)
dev.off()
wilcox.test(replicated.risk.genes.OSMR.results.summary$ribo.protein.slope[both.pass.index], replicated.risk.genes.OSMR.results.summary$ribo.protein.slope[up.pass.index])
##
png('se.downstream.png')
boxplot(replicated.risk.genes.OSMR.results.summary$se.of.ribo.protein.slope[both.pass.index], replicated.risk.genes.OSMR.results.summary$se.of.ribo.protein.slope[up.pass.index], names=c("pass", "fail"), main="Downstream MR: se of slope", outline = F)
dev.off()
wilcox.test(replicated.risk.genes.OSMR.results.summary$se.of.ribo.protein.slope[both.pass.index], replicated.risk.genes.OSMR.results.summary$se.of.ribo.protein.slope[up.pass.index])
##
png('f-stat.downstream.png')
boxplot(replicated.risk.genes.OSMR.results.summary$statistic.ribo[both.pass.index], replicated.risk.genes.OSMR.results.summary$statistic.ribo[up.pass.index], names=c("pass", "fail"), main="Downstream MR: instrument F statistics", outline = F)
dev.off()
wilcox.test(replicated.risk.genes.OSMR.results.summary$statistic.ribo[both.pass.index], replicated.risk.genes.OSMR.results.summary$statistic.ribo[up.pass.index])

## comparing both RNA->ribo and ribo -> protein MR stats  between both pass gene and no pass genes
both.pass.slope <- c(replicated.risk.genes.OSMR.results.summary$ribo.protein.slope[both.pass.index],replicated.risk.genes.OSMR.results.summary$rna.ribo.slope[both.pass.index])
no.pass.slope <- c(replicated.risk.genes.OSMR.results.summary$ribo.protein.slope[no.pass.index],replicated.risk.genes.OSMR.results.summary$rna.ribo.slope[no.pass.index])
##
png('effect_size.mr.png')
boxplot(both.pass.slope, no.pass.slope, names = c("pass", "fail"), main="MR: causal effect slope", outline = F)
dev.off()
wilcox.test(both.pass.slope, no.pass.slope)
##
both.pass.se <- c(replicated.risk.genes.OSMR.results.summary$se.of.ribo.protein.slope[both.pass.index],replicated.risk.genes.OSMR.results.summary$se.of.rna.ribo.slope[both.pass.index]) 
no.pass.se <- c(replicated.risk.genes.OSMR.results.summary$se.of.ribo.protein.slope[no.pass.index],replicated.risk.genes.OSMR.results.summary$se.of.rna.ribo.slope[no.pass.index])
png('se.mr.png')
boxplot(both.pass.se, no.pass.se, names = c("pass", "fail"), main="MR: se of slope", outline = F)
dev.off()
wilcox.test(both.pass.se, no.pass.se)
##
both.pass.F <- c(replicated.risk.genes.OSMR.results.summary$statistic.ribo[both.pass.index],replicated.risk.genes.OSMR.results.summary$statistic.RNA[both.pass.index]) 
no.pass.F <- c(replicated.risk.genes.OSMR.results.summary$statistic.ribo[no.pass.index],replicated.risk.genes.OSMR.results.summary$statistic.RNA[no.pass.index])
png('f-stat.mr.png')
boxplot(both.pass.F, no.pass.F, names = c("pass", "fail"), main="MR: instrument F statistics", outline = F)
dev.off()
wilcox.test(both.pass.F, no.pass.F)


