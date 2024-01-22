.libPaths("/zs32/home/tychu/R/lib64/R/library")
Args <- commandArgs()
#library(RColorBrewer)
#display.brewer.all(type = "qual")
#colors = brewer.pal(9,"Pastel2")[2:4]
#colors
#library(CMplot)
library(dplyr);library(qqman)
library( purrr );library(tidyr)
library(ggplot2);library(ggrepel)
library(optparse)
rna<-read.table(paste0("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/both_pathway/ENSG00000114904/eqtl.ENSG00000114904.txt"),sep="\t",header=T,comment.char="")
ribo<-read.table(paste0("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/both_pathway/ENSG00000114904/rqtl.ENSG00000114904.txt"),sep="\t",header=T,comment.char="")
protein<-read.table(paste0("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/both_pathway/ENSG00000114904/pqtl.ENSG00000114904.txt"),sep="\t",header=T,comment.char="")
gwas<-read.table(paste0("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/both_pathway/ENSG00000114904/gwas.ENSG00000114904.txt"),sep="\t",header=T,comment.char="")
head(rna)
rna<-rna[,c('X.chrom','pos','snp','pvalue')];colnames(rna)<-c('chr','position','snp','pvalue')
ribo<-ribo[,c('X.chrom','pos','snp','pvalue')];colnames(ribo)<-c('chr','position','snp','pvalue')
protein<-protein[,c('X.chrom','pos','snp','pvalue')];colnames(protein)<-c('chr','position','snp','pvalue')
gwas<-gwas[,c('X.chrom','pos','snp','pvalue')];colnames(gwas)<-c('chr','position','snp','pvalue')
## top snp
#coloc_top_snp_pos<-rna[order(rna$pvalue,decreasing=F),][1,]$position
## top coloc PP4 snp
coloc_top_snp_data<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/both_pathway/ENSG00000114904/top_coloc_snp')
gene<-'ENSG00000114904'
name<-'NEK4'
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
   geom_point(alpha=0.8, size=1.3,shape=1,color="lightblue4")+ 
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   #coord_cartesian(clip = 'off',ylim = c(0,15))+
   #annotate("text",x=52500000,y=22,label=coloc_top_snp,size=8,col="red")+
   ylim(0,pos_max)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
#   geom_label_repel( data=subset(Snp_pos_rna, gene_annotate=="yes"), aes(label=genename),size=4, col = "black",label.size=NA,force=1,fill = alpha(c("white"),0.9)) +
 #  geom_label_repel(data=all_data_rna,aes(label=genename),box.padding=unit(0.35, "lines"),size=4, col = "white",fill="red",fontface="bold",nudge_y=0.5,nudge_x=1.8)+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
   panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_ribo <- ggplot(ribo, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="lightblue4") + 
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   ylim(0,pos_max)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_protein <- ggplot(protein, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="lightblue4")+
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+ 
   ylim(0,pos_max)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_gwas <- ggplot(gwas, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="lightblue4")+
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   ylim(0,pos_max)+
   theme_bw() + labs(x='Position',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
library(cowplot)
library(patchwork)
library(ggplotify)
#p_mah<-plot_grid(title,p_mah_rna,p_mah_ribo,p_mah_protein,p_mah_gwas,nrow=5,align = 'v', axis = 'l',hjust=c(0,-12,-3,-10.5,-12),labels=c('','mRNA','ribosome occupancy','protein','GWAS'),label_size=12)
p_mah<-as.ggplot(plot_grid(p_mah_rna,p_mah_ribo,p_mah_protein,p_mah_gwas,nrow=4,align = 'v', hjust=c(-11.6,-3,-10.3,-11.6),labels=c('mRNA','ribosome occupancy','protein','GWAS'),label_size=16)+ plot_annotation(title = coloc_top_snp,theme = theme(plot.title = element_text(size = 22,hjust=0.5,color='red'))))
######### boxplot
library(ggplot2)
library(ggpubr)
data<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/both_pathway/ENSG00000114904/all.txt'),header=T,sep='\t')
data$genotype<-as.factor(data$genotype)
data$qtl<-factor(data$qtl,levels=c('mRNA','Ribo','Protein'),labels=c('mRNA','Ribo','Protein'))
head(data)
qtl_top_snp_data<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/both_pathway/ENSG00000114904/top_eqtl_snp')
qtl_top_snp=qtl_top_snp_data[1,1]
p_box<-ggboxplot(data, x = "qtl", y = "expression",
          color = "genotype",add = "jitter", palette = "jama")+
   guides(color = guide_legend(title = paste0(qtl_top_snp,' genotype')))+
   theme(legend.title=element_text(size=24),legend.text=element_text(size=20),plot.margin = unit(c(1,0,4,0),"lines"))+
   labs(x=NULL,title='NEK4')+theme(plot.title=element_text(size=30,face="italic"),axis.title.x=element_text(size=24),axis.text.x=element_text(size=20),axis.title.y=element_text(size=24),axis.text.y=element_text(size=20))
#design<-'#AA#BB#
         #AA#BB#
         #AA#BB#'
#p1<-wrap_plots(A=p_box,B=p_mah,design=design,widths=unit(c(1,10,10,1,14,14,2),c('cm','cm','cm','cm','cm','cm','cm')),heights=unit(c(12,7,1),c('cm','cm','cm')))
#+plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20,face='bold'))
#p1<-p1+plot_annotation(title = "CCDC117") &  theme(plot.title = element_text(hjust = 0.5,size=30,face="italic"))

###### no pathway
#no_gene<-Args[6]
#no_name<-Args[7]
no_gene<-"ENSG00000115524"
no_name<-"SF3B1"
rna<-read.table(paste0("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/no_pathway/",no_gene,"/eqtl.",no_gene,".txt"),sep="\t",header=T,comment.char="")
ribo<-read.table(paste0("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/no_pathway/",no_gene,"/rqtl.",no_gene,".txt"),sep="\t",header=T,comment.char="")
protein<-read.table(paste0("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/no_pathway/",no_gene,"/pqtl.",no_gene,".txt"),sep="\t",header=T,comment.char="")
gwas<-read.table(paste0("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/no_pathway/",no_gene,"/gwas.",no_gene,".txt"),sep="\t",header=T,comment.char="")
head(rna)
rna<-rna[,c('X.chrom','pos','snp','pvalue')];colnames(rna)<-c('chr','position','snp','pvalue')
ribo<-ribo[,c('X.chrom','pos','snp','pvalue')];colnames(ribo)<-c('chr','position','snp','pvalue')
protein<-protein[,c('X.chrom','pos','snp','pvalue')];colnames(protein)<-c('chr','position','snp','pvalue')
gwas<-gwas[,c('X.chrom','pos','snp','pvalue')];colnames(gwas)<-c('chr','position','snp','pvalue')
## top snp
#coloc_top_snp_pos<-rna[order(rna$pvalue,decreasing=F),][1,]$position
## top coloc PP4 snp
coloc_top_snp_data<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/no_pathway/',no_gene,'/top_coloc_snp'))
coloc_top_snp<-coloc_top_snp_data[1,1]
coloc_top_snp_pos<-as.numeric(strsplit(coloc_top_snp,"_")[[1]][2])
coloc_top_snp_pos
rna_pos_max<-max(-log10(rna$pvalue))
ribo_pos_max<-max(-log10(ribo$pvalue))
protein_pos_max<-max(-log10(protein$pvalue))
gwas_pos_max<-max(-log10(gwas$pvalue))
pos_max<-ceiling(max(c(rna_pos_max,ribo_pos_max,protein_pos_max,gwas_pos_max)))
pos_max

p_mah_rna <- ggplot(rna, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="lightblue4")+
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   #coord_cartesian(clip = 'off',ylim = c(0,6))+
   #annotate("text",x=40200000,y=8,label=coloc_top_snp,size=8,col="red")+
   ylim(0,pos_max)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
#   geom_label_repel( data=subset(Snp_pos_rna, gene_annotate=="yes"), aes(label=genename),size=4, col = "black",label.size=NA,force=1,fill = alpha(c("white"),0.9)) +
 #  geom_label_repel(data=all_data_rna,aes(label=genename),box.padding=unit(0.35, "lines"),size=4, col = "white",fill="red",fontface="bold",nudge_y=0.5,nudge_x=1.8)+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
   panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_ribo <- ggplot(ribo, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="lightblue4") +
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   ylim(0,pos_max)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_protein <- ggplot(protein, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="lightblue4")+
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   ylim(0,pos_max)+
   theme_bw() + labs(x='',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah_gwas <- ggplot(gwas, aes(x=position, y=-log10(pvalue))) +
   geom_point(alpha=0.8, size=1.3,shape=1,color="lightblue4")+
   geom_vline(aes(xintercept=coloc_top_snp_pos),color='red')+
   ylim(0,pos_max)+
   theme_bw() + labs(x='Position',y=expression(-log[10]*(p)))+
   theme(axis.title.x=element_text(size=20),axis.text.x=element_text(size=18),axis.title.y=element_text(size=20),axis.text.y=element_text(size=18))+
   theme(legend.position="none",panel.border = element_blank(),axis.line.y = element_line(),plot.margin = unit(c(0,2,0,0),"lines"),
     panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p_mah2<-as.ggplot(plot_grid(p_mah_rna,p_mah_ribo,p_mah_protein,p_mah_gwas,nrow=4,align = 'v', hjust=c(-11.6,-3,-10.3,-11.6),labels=c('mRNA','ribosome occupancy','protein','GWAS'),label_size=16)+ plot_annotation(title = coloc_top_snp,theme = theme(plot.title = element_text(size = 22,hjust=0.5,color='red'))))
######### boxplot
library(ggplot2)
library(ggpubr)
data<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/no_pathway/',no_gene,'/all.txt'),header=T,sep='\t')
data$genotype<-as.factor(data$genotype)
data$qtl<-factor(data$qtl,levels=c('mRNA','Ribo','Protein'),labels=c('mRNA','Ribo','Protein'))
head(data)
qtl_top_snp_data<-read.table(paste0('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/4.mr/fig4/no_pathway/',no_gene,'/top_eqtl_snp'))
qtl_top_snp=qtl_top_snp_data[1,1]
qtl_top_snp
p_box2<-ggboxplot(data, x = "qtl", y = "expression",
        color = "genotype",add = "jitter", palette = "jama")+
        guides(color = guide_legend(title = paste0(qtl_top_snp,' genotype')))+
   theme(legend.title=element_text(size=24),legend.position="top",legend.text=element_text(size=20),plot.margin = unit(c(1,0,4,0),"lines"))+
   labs(x=NULL,title=no_name)+theme(plot.title=element_text(size=30,face="italic"),axis.title.x=element_text(size=24),axis.text.x=element_text(size=20),axis.title.y=element_text(size=24),axis.text.y=element_text(size=20))
#p2<-wrap_plots(A=p_box2,B=p_mah2,design=design,widths=unit(c(10,10,2,14,14,2),c('cm','cm','cm','cm','cm','cm')),heights=unit(c(12,7,1),c('cm','cm','cm')))
#+plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20,face='bold'))
#p2<-p2+plot_annotation(title = "RAB15") &  theme(plot.title = element_text(hjust = 0.5,size=30,face="italic"))
#########
design<-'#AA#BB#
         #AA#BB#
         #CC#DD#
         #CC#DD#'
p<-wrap_plots(A=p_box,B=p_mah,C=p_box2,D=p_mah2,design=design,widths=unit(c(2,11,10,2,12,12,2),c('cm','cm','cm','cm','cm','cm','cm')),heights=unit(c(12,4,12,4),c('cm','cm','cm','cm')))+plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 30,face='bold'))
#p<-plot_grid(p_box,p_ ncol = 1, align = 'v',labels = c('A','B'),label_size = 30)
ggsave(paste0('fig4.png'),p,width=22,height=18,bg='white')

