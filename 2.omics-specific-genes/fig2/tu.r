.libPaths("/zs32/home/tychu/R/lib64/R/library")
Args <- commandArgs()
library(dplyr);library(qqman)
library( purrr );library(tidyr)
library(ggplot2);library(ggrepel)
library(optparse)
rna<-read.table(paste0("esQTL/ENSG00000159873/eqtl.ENSG00000159873.txt"),sep="\t",header=T,comment.char="")
ribo<-read.table(paste0("esQTL/ENSG00000159873/rqtl.ENSG00000159873.txt"),sep="\t",header=T,comment.char="")
protein<-read.table(paste0("esQTL/ENSG00000159873/pqtl.ENSG00000159873.txt"),sep="\t",header=T,comment.char="")
gwas<-read.table(paste0("esQTL/ENSG00000159873/gwas.ENSG00000159873.txt"),sep="\t",header=T,comment.char="")
head(rna)
rna<-rna[,c('X.chrom','pos','snp','pvalue')];colnames(rna)<-c('chr','position','snp','pvalue')
ribo<-ribo[,c('X.chrom','pos','snp','pvalue')];colnames(ribo)<-c('chr','position','snp','pvalue')
protein<-protein[,c('X.chrom','pos','snp','pvalue')];colnames(protein)<-c('chr','position','snp','pvalue')
gwas<-gwas[,c('X.chrom','pos','snp','pvalue')];colnames(gwas)<-c('chr','position','snp','pvalue')
## top snp
coloc_top_snp_data<-read.table('esQTL/ENSG00000159873/top_coloc_snp')
coloc_top_snp<-coloc_top_snp_data[1,1]
coloc_top_snp
coloc_top_snp_pos<-as.numeric(strsplit(coloc_top_snp,"_")[[1]][2])
coloc_top_snp_pos
rna_pos_max<-max(-log10(rna$pvalue))
ribo_pos_max<-max(-log10(ribo$pvalue))
protein_pos_max<-max(-log10(protein$pvalue))
gwas_pos_max<-max(-log10(gwas$pvalue))
pos_max<-ceiling(max(c(rna_pos_max,ribo_pos_max,protein_pos_max,gwas_pos_max)))
pos_max

##
p_mah_rna <- ggplot(rna, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="azure4")+ 
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   ylim(0,6)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(plot.title=element_text(size=20,hjust=0.5,color='red'))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
   panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_ribo <- ggplot(ribo, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="azure4") + 
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   ylim(0,6)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_protein <- ggplot(protein, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="azure4")+
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+ 
   ylim(0,6)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_gwas <- ggplot(gwas, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="azure4")+
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   ylim(0,6)+
   theme_bw() + labs(x='Position',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
library(cowplot)
library(patchwork)
library(ggplotify)
p_mah<-as.ggplot(plot_grid(p_mah_rna,p_mah_ribo,p_mah_protein,p_mah_gwas,nrow=4,align = 'v',hjust=c(-11.5,-3.2,-10.3,-11.5),labels=c('mRNA','ribosome occupancy','protein','GWAS'),label_size=16) + plot_annotation(title = coloc_top_snp,theme = theme(plot.title = element_text(size = 20,hjust=0.5,color='red'))))
######### boxplot
library(ggplot2)
library(ggpubr)
data<-read.table(paste0('esQTL/ENSG00000159873/all.txt'),header=T,sep='\t')
data$genotype<-as.factor(data$genotype)
data$qtl<-factor(data$qtl,levels=c('mRNA','Ribo','Protein'),labels=c('mRNA','ribosome\noccupancy','protein'))
head(data)
qtl_top_snp_data<-read.table(paste0('esQTL/ENSG00000159873/top_eqtl_snp'))
qtl_top_snp=qtl_top_snp_data[1,1]
qtl_top_snp
p_box<-ggboxplot(data, x = "qtl", y = "expression",
          color = "genotype",add = "jitter", palette = "jama")+
   guides(color = guide_legend(title = paste0(qtl_top_snp,' genotype')))+
   theme(legend.title=element_text(size=24),legend.text=element_text(size=20),plot.margin = unit(c(1,0,4,0),"lines"))+
   labs(x=NULL,title='CCDC117')+theme(plot.title=element_text(size=30,face="italic"),axis.title.x=element_text(size=24),axis.text.x=element_text(size=20),axis.title.y=element_text(size=24),axis.text.y=element_text(size=20))

######
rsqtl_gene<-Args[6]
rsqtl_name<-Args[7]
rna<-read.table(paste0('rsQTL/',rsqtl_gene,'/eqtl.',rsqtl_gene,'.txt'),sep="\t",header=T,comment.char="")
ribo<-read.table(paste0('rsQTL/',rsqtl_gene,'/rqtl.',rsqtl_gene,'.txt'),sep="\t",header=T,comment.char="")
protein<-read.table(paste0('rsQTL/',rsqtl_gene,'/pqtl.',rsqtl_gene,'.txt'),sep="\t",header=T,comment.char="")
gwas<-read.table(paste0('rsQTL/',rsqtl_gene,'/gwas.',rsqtl_gene,'.txt'),sep="\t",header=T,comment.char="")
head(rna)
rna<-rna[,c('X.chrom','pos','snp','pvalue')];colnames(rna)<-c('chr','position','snp','pvalue')
ribo<-ribo[,c('X.chrom','pos','snp','pvalue')];colnames(ribo)<-c('chr','position','snp','pvalue')
protein<-protein[,c('X.chrom','pos','snp','pvalue')];colnames(protein)<-c('chr','position','snp','pvalue')
gwas<-gwas[,c('X.chrom','pos','snp','pvalue')];colnames(gwas)<-c('chr','position','snp','pvalue')
## top coloc PP4 snp
coloc_top_snp_data<-read.table(paste0('rsQTL/',rsqtl_gene,'/top_coloc_snp'))
coloc_top_snp<-coloc_top_snp_data[1,1]
coloc_top_snp_pos<-as.numeric(strsplit(coloc_top_snp,"_")[[1]][2])
coloc_top_snp_pos
coloc_top_snp
rna_pos_max<-max(-log10(rna$pvalue))
ribo_pos_max<-max(-log10(ribo$pvalue))
protein_pos_max<-max(-log10(protein$pvalue))
gwas_pos_max<-max(-log10(gwas$pvalue))
pos_max<-ceiling(max(c(rna_pos_max,ribo_pos_max,protein_pos_max,gwas_pos_max)))
pos_max

p_mah_rna <- ggplot(rna, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="azure4")+
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   ylim(0,20)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(plot.title=element_text(size=20,hjust=0.5,color='red'))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
   panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_ribo <- ggplot(ribo, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="azure4") +
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_protein <- ggplot(protein, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="azure4")+
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   ylim(0,6)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_gwas <- ggplot(gwas, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="azure4")+
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   theme_bw() + labs(x='Position',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah2<-as.ggplot(plot_grid(p_mah_rna,p_mah_ribo,p_mah_protein,p_mah_gwas,nrow=4,align = 'v',hjust=c(-11.5,-3.2,-10.3,-11.5),labels=c('mRNA','ribosome occupancy','protein','GWAS'),label_size=16)+ plot_annotation(title = coloc_top_snp,theme = theme(plot.title = element_text(size = 20,hjust=0.5,color='red'))))
######### boxplot
library(ggplot2)
library(ggpubr)
data<-read.table(paste0('rsQTL/',rsqtl_gene,'/all.txt'),header=T,sep='\t')
data$genotype<-as.factor(data$genotype)
data$qtl<-factor(data$qtl,levels=c('mRNA','Ribo','Protein'),labels=c('mRNA','ribosome\noccupancy','protein'))
head(data)
qtl_top_snp_data<-read.table(paste0('rsQTL/',rsqtl_gene,'/top_rqtl_snp'))
qtl_top_snp=qtl_top_snp_data[1,1]
qtl_top_snp
p_box2<-ggboxplot(data, x = "qtl", y = "expression",
        color = "genotype",add = "jitter", palette = "jama")+
        guides(color = guide_legend(title = paste0(qtl_top_snp,' genotype')))+
   theme(legend.title=element_text(size=24),legend.position="top",legend.text=element_text(size=20),plot.margin = unit(c(1,0,4,0),"lines"))+
   labs(x=NULL,title=rsqtl_name)+theme(plot.title=element_text(size=30,face="italic"),axis.title.x=element_text(size=24),axis.text.x=element_text(size=20),axis.title.y=element_text(size=24),axis.text.y=element_text(size=20))
#########
design<-'#AA#BB#
         #AA#BB#
         #CC#DD#
         #CC#DD#'
p<-wrap_plots(A=p_box,B=p_mah,C=p_box2,D=p_mah2,design=design,widths=unit(c(1,11,10,2,14,12,1),c('cm','cm','cm','cm','cm','cm','cm')),heights=unit(c(12,4,12,4),c('cm','cm','cm','cm')))+plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 30,face='bold'))
ggsave(paste0('fig2.png'),p,width=21,height=19,bg='white')

