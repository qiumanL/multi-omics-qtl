.libPaths("/zs32/home/tychu/R/lib64/R/library")
library(VennDiagram)
library(ggVennDiagram)
library(ggplot2)
library(Hmisc)
Args <- commandArgs()
rna<-read.table('../rna/2021PGC3.asso.csv.coding.multiTest.bonferroni005.nonMHC.genename')
ribo<-read.table('../ribo/2021PGC3.asso.csv.coding.multiTest.bonferroni005.nonMHC.genename')
protein<-read.table('../protein/2021PGC3.asso.csv.gene.omnibus.multiTest.bonferroni005.nonMHC.genename')
library(RColorBrewer)
display.brewer.all(type = "qual")
colors = brewer.pal(9,"Pastel2")[2:4]
colors
x=list(mRNA=rna[,1],Ribo_seq=ribo[,1],Mass_spec=protein[,1])
#png('venn.risk.png')
p_venn<-ggVennDiagram(x,label="count",label_percent_digit=2,label_size=8,label_alpha = 0,
  category.names = c("mRNA","ribosome\noccupancy","protein"), set_size = 5.5)+
  scale_fill_gradient(low="white",high = "lightblue")+
  scale_color_manual(values=colors)+
  scale_x_continuous(expand = expansion(mult = .1))+
  guides(fill=F)
#dev.off()
print("detected genes at 3 omics:")
intersect(rna[,2],intersect(ribo[,2],protein[,2]))
print("detected genes at only rna and ribo:")
setdiff(intersect(rna[,2],ribo[,2]),intersect(rna[,2],intersect(ribo[,2],protein[,2])))
print("detected genes at only rna and protein:")
setdiff(intersect(rna[,2],protein[,2]),intersect(rna[,2],intersect(ribo[,2],protein[,2])))
print("detected genes at only ribo and protein:")
setdiff(intersect(ribo[,2],protein[,2]),intersect(rna[,2],intersect(ribo[,2],protein[,2])))
print("detected genes at only rna:")
setdiff(rna[,2],union(ribo[,2],protein[,2]))
print("detected genes at only ribo:")
setdiff(ribo[,2],union(rna[,2],protein[,2]))
print("detected genes at only protein:")
setdiff(protein[,2],union(rna[,2],ribo[,2]))
##################
library(dplyr);library(qqman)
library( purrr );library(tidyr)
library(ggplot2);library(ggrepel)
library(optparse)
rna<-read.table("../rna/2021PGC3.asso.csv.coding.multiTest.mahattan",sep="\t",header=F)
ribo<-read.table("../ribo/2021PGC3.asso.csv.coding.multiTest.mahattan",sep="\t",header=F)
protein<-read.table("../protein/2021PGC3.asso.csv.gene.omnibus.multiTest.mahattan",sep="\t",header=F)
colnames(rna)<-c('gene','chr','position','rna','genename')
colnames(ribo)<-c('gene','chr','position','ribo','genename')
colnames(protein)<-c('gene','chr','position','protein','genename')
risk_gene_rna<-read.table('../rna/2021PGC3.asso.csv.coding.multiTest.bonferroni005.nonMHC',sep="\t",header=F)[,1]
risk_gene_ribo<-read.table('../ribo/2021PGC3.asso.csv.coding.multiTest.bonferroni005.nonMHC',sep="\t",header=F)[,1]
risk_gene_protein<-read.table('../protein/2021PGC3.asso.csv.gene.omnibus.multiTest.bonferroni005.nonMHC',sep="\t",header=F)[,1]
rna<- rna %>% mutate(gene_annotate=ifelse(gene %in% risk_gene_rna,"yes","no"))
rna<-subset(rna,!(rna<0.05 & gene_annotate=='no'))
chr_len_rna <- rna %>% group_by(chr) %>% summarise(chr_len=max(position))
head(chr_len_rna)
chr_pos_rna <- chr_len_rna  %>% mutate(total = cumsum(chr_len) - chr_len) %>% select(-chr_len)
head(chr_pos_rna)
Snp_pos_rna <- chr_pos_rna %>% left_join(rna, ., by="chr") %>% arrange(chr, position) %>% mutate( BPcum = position + total)
X_axis_rna <-  Snp_pos_rna %>% group_by(chr) %>% summarize(center=( max(BPcum) +min(BPcum) ) / 2 )
######
ribo <- ribo %>% mutate(gene_annotate=ifelse(gene %in% risk_gene_ribo,"yes","no"))
ribo<-subset(ribo,!(ribo<0.05 & gene_annotate=='no'))
chr_len_ribo <- ribo %>% group_by(chr) %>% summarise(chr_len=max(position))
chr_pos_ribo <- chr_len_ribo  %>% mutate(total = cumsum(chr_len) - chr_len) %>% select(-chr_len)
Snp_pos_ribo <- chr_pos_ribo %>% left_join(ribo, ., by="chr") %>% arrange(chr, position) %>% mutate( BPcum = position + total)
X_axis_ribo <-  Snp_pos_ribo %>% group_by(chr) %>% summarize(center=( max(BPcum) +min(BPcum) ) / 2 )
head(Snp_pos_ribo)
######
protein <- protein %>% mutate(gene_annotate=ifelse(gene %in% risk_gene_protein,"yes","no"))
protein<-subset(protein,!(protein<0.05 & gene_annotate=='no'))
chr_len_protein <- protein %>% group_by(chr) %>% summarise(chr_len=max(position))
chr_pos_protein <- chr_len_protein  %>% mutate(total = cumsum(chr_len) - chr_len) %>% select(-chr_len)
Snp_pos_protein <- chr_pos_protein %>% left_join(protein, ., by="chr") %>% arrange(chr, position) %>% mutate( BPcum = position + total)
head(Snp_pos_protein)
X_axis_protein <-  Snp_pos_protein %>% group_by(chr) %>% summarize(center=( max(BPcum) +min(BPcum) ) / 2 )
###
both<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/replicate_summary/both_pathways.txt',header=T,sep='\t')
one<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/replicate_summary/rna_ribo_pathway.txt',header=T,sep='\t')
no<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/replicate_summary/no_pathway.txt',header=T,sep='\t')
Snp_pos_rna_risk<-subset(Snp_pos_rna,gene_annotate=="yes")
Snp_pos_ribo_risk<-subset(Snp_pos_ribo,gene_annotate=="yes")
Snp_pos_protein_risk<-subset(Snp_pos_protein,gene_annotate=="yes")
print("risk gene num")
dim(Snp_pos_rna_risk)
dim(Snp_pos_ribo_risk)
dim(Snp_pos_protein_risk)
## no risk
norisk_gene_rna<-subset(Snp_pos_rna,gene_annotate=="no")$genename
norisk_gene_ribo<-subset(Snp_pos_ribo,gene_annotate=="no")$genename
norisk_gene_protein<-subset(Snp_pos_protein,gene_annotate=="no")$genename

## replicate
rep_gene_rna<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/replicate_summary/scz.eQTL.MR.results.rep.txt',header=T,sep='\t')$gene.name
rep_gene_ribo<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/replicate_summary/scz.rQTL.MR.results.rep.txt',header=T,sep='\t')$gene.name
rep_gene_protein<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/replicate_summary/scz.pQTL.MR.results.rep.txt',header=T,sep='\t')$gene.name
##
Snp_pos_rna_rep<-Snp_pos_rna[Snp_pos_rna$genename %in% rep_gene_rna,]
Snp_pos_ribo_rep<-Snp_pos_ribo[Snp_pos_ribo$genename %in% rep_gene_ribo,]
Snp_pos_protein_rep<-Snp_pos_protein[Snp_pos_protein$genename %in% rep_gene_protein,]
print("replicated risk gene num")
dim(Snp_pos_rna_rep)
dim(Snp_pos_ribo_rep)
dim(Snp_pos_protein_rep)
print("#######")
##
norep_data_rna<-Snp_pos_rna_risk[-which(Snp_pos_rna_risk$genename %in% rep_gene_rna),]
norep_data_ribo<-Snp_pos_ribo_risk[-which(Snp_pos_ribo_risk$genename %in% rep_gene_ribo),]
norep_data_protein<-Snp_pos_protein_risk[-which(Snp_pos_protein_risk$genename %in% rep_gene_protein),]
norep_gene_rna<-norep_data_rna$genename
norep_gene_ribo<-norep_data_ribo$genename
norep_gene_protein<-norep_data_protein$genename
print("not replicated risk gene num")
dim(norep_data_rna)
dim(norep_data_ribo)
dim(norep_data_protein)
print("#######")
##
both_gene<-both$gene.name.RNAtoRIBO
one_gene<-one$gene.name.RNAtoRIBO
no_gene<-no$gene.name.RNAtoRIBO
##
Snp_pos_rna$type<-'control'
Snp_pos_rna$type[which(Snp_pos_rna$genename %in% norisk_gene_rna & Snp_pos_rna$chr %% 2 == 1)]<-'norisk_odd'
Snp_pos_rna$type[which(Snp_pos_rna$genename %in% norisk_gene_rna & Snp_pos_rna$chr %% 2 == 0)]<-'norisk_even'
Snp_pos_rna$type[which(Snp_pos_rna$genename %in% norep_gene_rna)]<-'norep'
Snp_pos_rna$type[which(Snp_pos_rna$genename %in% both_gene & Snp_pos_rna$type %nin% c('norisk_odd','norisk_even','norep'))]<-'both'
Snp_pos_rna$type[which(Snp_pos_rna$genename %in% one_gene & Snp_pos_rna$type %nin% c('norisk_odd','norisk_even','norep'))]<-'one'
Snp_pos_rna$type[which(Snp_pos_rna$genename %in% no_gene & Snp_pos_rna$type %nin% c('norisk_odd','norisk_even','norep'))]<-'no'
table(Snp_pos_rna$type)
Snp_pos_rna$type<-factor(Snp_pos_rna$type,levels=c('norisk_odd','norisk_even','norep','both','one','no'),labels=c('norisk_odd','norisk_even','norep','both','one','no'))
Snp_pos_rna_risk<-Snp_pos_rna[Snp_pos_rna$type %in% c('norep','both','one','no'),]
table(Snp_pos_rna_risk$type)
##
Snp_pos_ribo$type<-'control'
Snp_pos_ribo$type[which(Snp_pos_ribo$genename %in% norisk_gene_ribo & Snp_pos_ribo$chr %% 2 == 1)]<-'norisk_odd'
Snp_pos_ribo$type[which(Snp_pos_ribo$genename %in% norisk_gene_ribo & Snp_pos_ribo$chr %% 2 == 0)]<-'norisk_even'
Snp_pos_ribo$type[which(Snp_pos_ribo$genename %in% norep_gene_ribo)]<-'norep'
Snp_pos_ribo$type[which(Snp_pos_ribo$genename %in% both_gene  & Snp_pos_ribo$type %nin% c('norisk_odd','norisk_even','norep'))]<-'both'
Snp_pos_ribo$type[which(Snp_pos_ribo$genename %in% one_gene & Snp_pos_ribo$type %nin% c('norisk_odd','norisk_even','norep'))]<-'one'
Snp_pos_ribo$type[which(Snp_pos_ribo$genename %in% no_gene & Snp_pos_ribo$type %nin% c('norisk_odd','norisk_even','norep'))]<-'no'
Snp_pos_ribo$type<-factor(Snp_pos_ribo$type,levels=c('norisk_odd','norisk_even','norep','both','one','no'),labels=c('norisk_odd','norisk_even','norep','both','one','no'))
dim(Snp_pos_ribo)
Snp_pos_ribo_risk<-Snp_pos_ribo[Snp_pos_ribo$type %in% c('norep','both','one','no'),]
table(Snp_pos_ribo$type)
Snp_pos_ribo[Snp_pos_ribo$type=='one',]
##
Snp_pos_protein$type<-'control'
Snp_pos_protein$type[which(Snp_pos_protein$genename %in% norisk_gene_protein & Snp_pos_protein$chr %% 2 == 1)]<-'norisk_odd'
Snp_pos_protein$type[which(Snp_pos_protein$genename %in% norisk_gene_protein & Snp_pos_protein$chr %% 2 == 0)]<-'norisk_even'
Snp_pos_protein$type[which(Snp_pos_protein$genename %in% norep_gene_protein)]<-'norep'
Snp_pos_protein$type[which(Snp_pos_protein$genename %in% both_gene & Snp_pos_protein$type %nin% c('norisk_odd','norisk_even','norep'))]<-'both'
Snp_pos_protein$type[which(Snp_pos_protein$genename %in% one_gene & Snp_pos_protein$type %nin% c('norisk_odd','norisk_even','norep'))]<-'one'
Snp_pos_protein$type[which(Snp_pos_protein$genename %in% no_gene & Snp_pos_protein$type %nin% c('norisk_odd','norisk_even','norep'))]<-'no'
Snp_pos_protein$type<-factor(Snp_pos_protein$type,levels=c('norisk_odd','norisk_even','norep','both','one','no'),labels=c('norisk_odd','norisk_even','norep','both','one','no'))
head(Snp_pos_protein)
Snp_pos_protein_risk<-Snp_pos_protein[Snp_pos_protein$type %in% c('norep','both','one','no'),]
table(Snp_pos_protein_risk$type)
colors_man<-c('darkgreen','purple','darkgrey','darkred','cornflowerblue','orange')
shapes_man<-c(1,1,8,16,17,15)
##
p_mah_rna <- ggplot() +
   geom_point(Snp_pos_rna, mapping=aes(x=BPcum, y=-log10(rna),color=factor(type),shape=type), alpha=0.8, size=2) +
#   scale_color_manual(values = rep(c("darkgreen", "purple"), 22 )) +
   scale_color_manual(values=colors_man)+
   scale_shape_manual(values =shapes_man)+
   scale_x_continuous(label = X_axis_rna$chr, breaks= X_axis_rna$center ) +
#   scale_y_continuous(expand = c(0, 0) ) + ylim(c(0,18))+ 
   geom_hline(yintercept = c(-log10(0.05)), color = c('black'),size = 1.2, linetype = c("twodash")) + 
   theme_bw() + labs(x='',y=expression(-log[10]*(adj.p)))+
   theme(axis.title.x=element_text(size=16),axis.text.x=element_text(size=12),axis.title.y=element_text(size=16),axis.text.y=element_text(size=12))+
   geom_label_repel( data=Snp_pos_rna_risk, aes(x=BPcum,y=-log10(rna),label=genename,color=factor(type)),size=6,label.size=NA,force=1,fill = alpha(c("white"),0.9),max.overlaps=30) +
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())

p_mah_ribo <- ggplot() +
   geom_point(Snp_pos_ribo, mapping=aes(x=BPcum, y=-log10(ribo),color=factor(type),shape=type), alpha=0.8, size=2) +
   scale_color_manual(values = colors_man) +
   scale_shape_manual(values =shapes_man)+
   scale_x_continuous(label = X_axis_ribo$chr, breaks= X_axis_ribo$center ) +
#   scale_y_continuous(expand = c(0, 0) ) + ylim(c(0,18))+
   geom_hline(yintercept = c(-log10(0.05)), color = c('black'),size = 1.2, linetype = c("twodash")) +
   theme_bw() + labs(x='',y=expression(-log[10]*(adj.p)))+
   theme(axis.title.x=element_text(size=16),axis.text.x=element_text(size=12),axis.title.y=element_text(size=16),axis.text.y=element_text(size=12))+
   geom_label_repel(data=Snp_pos_ribo_risk, aes(x=BPcum,y=-log10(ribo),label=genename,color=factor(type)),size=6,label.size=NA,force=1,fill = alpha(c("white"),0.9),max.overlaps=30) +
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
print('ye')
####
p_mah_protein <- ggplot() +
   geom_point(Snp_pos_protein, mapping=aes(x=BPcum, y=-log10(protein),color=factor(type),shape=type), alpha=0.8, size=2) +
#   scale_color_manual(values = rep(c("darkgreen", "purple"), 22 )) +
   scale_color_manual(values =colors_man)+
   scale_shape_manual(values =shapes_man)+
   scale_x_continuous(label = X_axis_protein$chr, breaks= X_axis_protein$center ) +
   scale_y_continuous(expand = c(0, 0) ) + ylim(c(0,18))+
   geom_hline(yintercept = c(-log10(0.05)), color = c('black'),size = 1.2, linetype = c("twodash")) +
#   geom_point(data=all_data_protein, color="red", size=1.6) +
#   geom_point(data=norep_data_protein, color="darkgrey",shape=8, size=2.6) +
#   geom_point(data=both_data_protein, color="darkred",shape=16, size=2.6) +
#   geom_point(data=one_data_protein, color="cornflowerblue",shape=17, size=2.6) +
#   geom_point(data=no_data_protein, color="orange",shape=15, size=2.6) +
   theme_bw() + labs(x='Chromosome',y=expression(-log[10]*(adj.p)))+
   theme(axis.title.x=element_text(size=16),axis.text.x=element_text(size=12),axis.title.y=element_text(size=16),axis.text.y=element_text(size=12))+
#   geom_label_repel( data=norep_data_protein, aes(label=genename),col='darkgrey',size=6,label.size=NA,force=1,fill = alpha(c("white"),0)) +
#   geom_label_repel( data=both_data_protein, aes(label=genename),col='darkred',size=6,label.size=NA,force=1,fill = alpha(c("white"),0),nudge_y=0.8-log10(no_data_protein$protein)) +
#   geom_label_repel( data=one_data_protein, aes(label=genename),col='cornflowerblue',size=6,label.size=NA,force=1,fill = alpha(c("white"),0)) +
#   geom_label_repel( data=no_data_protein, aes(label=genename),col='orange',size=6,label.size=NA,force=1,fill = alpha(c("white"),0),nudge_y=0.5-log10(no_data_protein$protein)) +
   geom_label_repel( data=Snp_pos_protein_risk, aes(x=BPcum,y=-log10(protein),label=genename,color=factor(type)),size=6,label.size=NA,force=1,fill = alpha(c("white"),0.9)) +
#   scale_color_manual(values = c('darkgrey','black','cornflowerblue','orange' )) +
#   geom_label_repel(data=subset(Snp_pos_protein, gene_annotate=="yes"), aes(label=genename),size=4, col = "black",label.size=NA,force=1,fill = alpha(c("white"),0.9)) +
#   geom_label_repel(data=all_data_protein,aes(label=genename),box.padding=unit(0.35, "lines"),fontface="bold",col="white",fill="red",nudge_y=0)+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
#########
library(cowplot)
p_mah_all<-plot_grid(p_mah_rna,p_mah_ribo,p_mah_protein,nrow=3,align = 'v', hjust=c(-20,-6,-18),vjust=c(-0.5,0,0),labels=c('mRNA','ribosome occupancy','protein'),label_size=20)
###########  fig 3C
.libPaths("/zs32/home/tychu/R/lib64/R/library")
coloc_top_snp_data<-read.table('example_gene/ENSG00000138030/top_coloc_snp')
Args <- commandArgs()
gene<-'ENSG00000138030'
name<-'example_gene'
coloc_top_snp<-coloc_top_snp_data[1,1]
coloc_top_snp
#library(RColorBrewer)
#display.brewer.all(type = "qual")
#colors = brewer.pal(9,"Pastel2")[2:4]
#colors
#library(CMplot)
library(dplyr);library(qqman)
library( purrr );library(tidyr)
library(ggplot2);library(ggrepel)
library(optparse)
rna<-read.table(paste0("example_gene/ENSG00000138030/eqtl.",gene,".txt"),sep="\t",header=T,comment.char="")
ribo<-read.table(paste0("example_gene/ENSG00000138030/rqtl.",gene,".txt"),sep="\t",header=T,comment.char="")
protein<-read.table(paste0("example_gene/ENSG00000138030/pqtl.",gene,".txt"),sep="\t",header=T,comment.char="")
gwas<-read.table(paste0("example_gene/ENSG00000138030/gwas.",gene,".txt"),sep="\t",header=T,comment.char="")
head(rna)
rna<-rna[,c('X.chrom','pos','snp','pvalue')];colnames(rna)<-c('chr','position','snp','pvalue')
ribo<-ribo[,c('X.chrom','pos','snp','pvalue')];colnames(ribo)<-c('chr','position','snp','pvalue')
protein<-protein[,c('X.chrom','pos','snp','pvalue')];colnames(protein)<-c('chr','position','snp','pvalue')
gwas<-gwas[,c('X.chrom','pos','snp','pvalue')];colnames(gwas)<-c('chr','position','snp','pvalue')
## top snp
#top_snp_pos<-rna[order(rna$pvalue,decreasing=F),][1,]$position
## top coloc PP4 snp
coloc_top_snp_pos<-as.numeric(strsplit(coloc_top_snp,"_")[[1]][2])
coloc_top_snp_pos
coloc_top_snp_pos
##
p_mah_rna <- ggplot(rna, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="azure4")+ 
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   #coord_cartesian(clip = 'off',ylim = c(0,6))+
   #annotate("text",x=27004882,y=8,label=coloc_top_snp,size=6,col="red")+
#   ylim(0,30)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=18),axis.text.x=element_text(size=15),axis.title.y=element_text(size=18),axis.text.y=element_text(size=15))+
#   geom_label_repel( data=subset(Snp_pos_rna, gene_annotate=="yes"), aes(label=genename),size=4, col = "black",label.size=NA,force=1,fill = alpha(c("white"),0.9)) +
 #  geom_label_repel(data=all_data_rna,aes(label=genename),box.padding=unit(0.35, "lines"),size=4, col = "white",fill="red",fontface="bold",nudge_y=0.5,nudge_x=1.8)+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
   panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_ribo <- ggplot(ribo, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="azure4") + 
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
 #  ylim(0,30)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=18),axis.text.x=element_text(size=15),axis.title.y=element_text(size=18),axis.text.y=element_text(size=15))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_protein <- ggplot(protein, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="azure4")+
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+ 
  # ylim(0,30)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=18),axis.text.x=element_text(size=15),axis.title.y=element_text(size=18),axis.text.y=element_text(size=15))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_gwas <- ggplot(gwas, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="azure4")+
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   #ylim(0,30)+
   theme_bw() + labs(x='Position',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=18),axis.text.x=element_text(size=15),axis.title.y=element_text(size=18),axis.text.y=element_text(size=15))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
library(cowplot)
library(patchwork)
library(ggplotify)
#p_mah<-plot_grid(title,p_mah_rna,p_mah_ribo,p_mah_protein,p_mah_gwas,nrow=5,align = 'v', axis = 'l',hjust=c(0,-12,-3,-10.5,-12),labels=c('','mRNA','ribosome occupancy','protein','GWAS'),label_size=12)
p_mah<-as.ggplot(plot_grid(p_mah_rna,p_mah_ribo,p_mah_protein,p_mah_gwas,nrow=4,align = 'v', hjust=c(-12,-3.2,-10.5,-12),vjust=c(0,0,0,-0.5),labels=c('mRNA','ribosome occupancy','protein','GWAS'),label_size=12)+ plot_annotation(title = coloc_top_snp,theme = theme(plot.title = element_text(size = 20,hjust=0.5,color='red'))))
######### boxplot
library(ggplot2)
library(ggpubr)
data<-read.table(paste0('example_gene/ENSG00000138030/all.txt'),header=T,sep='\t')
qtl_top_snp_data<-read.table('example_gene/ENSG00000138030/top_eqtl_snp')
data$genotype<-as.factor(data$genotype)
data$qtl<-factor(data$qtl,levels=c('mRNA','Ribo','Protein'),labels=c('mRNA','ribosome\noccupancy','protein'))
head(data)
qtl_top_snp=qtl_top_snp_data[1,1]
qtl_top_snp
p_box<-ggboxplot(data, x = "qtl", y = "expression",
   color = "genotype",add = "jitter", palette = "jama")+
   guides(color = guide_legend(title = paste0(qtl_top_snp,' genotype')))+
   theme(legend.position = "top")+
   theme(legend.title=element_text(size=18),legend.text=element_text(size=15),plot.margin = unit(c(1,0,4,0),"lines"))+
   labs(x=NULL,title=name)+theme(axis.title.x=element_text(size=18),axis.text.x=element_text(size=15),axis.title.y=element_text(size=18),axis.text.y=element_text(size=15))+
   theme(plot.title=element_text(face="italic",size=22))
p<-plot_grid(p_box,p_mah)
design<-'##########
         #AAAAAAAA#
         #AAAAAAAA#
         ##########
         #BB#CC#DD#
         #BB#CC#DD#
         #BB#CC#DD#
         #BB#CC#DD#'
p2<-wrap_plots(A=p_mah_all,B=p_venn,C=p_box,D=p_mah,design=design,widths=unit(c(1,6,8,1,6,8,1,10,10,2),c('cm','cm','cm','cm','cm','cm','cm','cm','cm','cm')),heights=unit(c(1,14,5,0.995,0.005,10,2,2),c('cm','cm','cm','cm','cm','cm','cm','cm')))+plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 27,face='bold'))
ggsave(paste0('Fig3.png'),p2,width=22,height=17,bg='white')

