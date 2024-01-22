library(ggplot2)
library(ggsci)
library(ggpattern)
library(ggplot2)
library(dplyr)
library(Rmisc)
eqtl<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/eQTL/eqtl.permute.qvalue.txt',sep="\t",head=T)
rqtl<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/rQTL/rqtl.permute.qvalue.txt',sep="\t",head=T)
pqtl<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/pQTL/uniq/pqtl.permute.qvalue.txt',sep="\t",head=T)
eqtl<-eqtl$P.beta
rqtl<-rqtl$P.beta
pqtl<-pqtl$P.beta
eqtl_sig<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/eQTL/eqtl.permute.qvalue.sig.txt',sep="\t",head=T)
rqtl_sig<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/rQTL/rqtl.permute.qvalue.sig.txt',sep="\t",head=T)
pqtl_sig<-read.table('/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/pQTL/pqtl.permute.qvalue.sig.txt',sep="\t",head=T)
library(tidyr)
eqtl_sig_num<-dim(eqtl_sig)[1]
rqtl_sig_num<-dim(rqtl_sig)[1]
pqtl_sig_num<-dim(pqtl_sig)[1]
eqtl_sig_num
rqtl_sig_num
pqtl_sig_num
n_eqtl  <- length(eqtl)
ci=0.95
#ppoints(n_eqtl)
df_eqtl <- data.frame(
    observed = -log10(sort(eqtl)),
    expected = -log10(ppoints(n_eqtl)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n_eqtl, shape2 = n_eqtl:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n_eqtl, shape2 = n_eqtl:1)),
    qtl = rep('eQTL',n_eqtl)
  )
n_rqtl  <- length(rqtl)
df_rqtl <- data.frame(
    observed = -log10(sort(rqtl)),
    expected = -log10(ppoints(n_rqtl)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n_rqtl, shape2 = n_rqtl:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n_rqtl, shape2 = n_rqtl:1)),
    qtl = rep('rQTL',n_rqtl)
  )
n_pqtl  <- length(pqtl)
df_pqtl <- data.frame(
    observed = -log10(sort(pqtl)),
    expected = -log10(ppoints(n_pqtl)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n_pqtl, shape2 = n_pqtl:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n_pqtl, shape2 = n_pqtl:1)),
    qtl = rep('pQTL',n_pqtl)
)
df<- rbind(df_eqtl,df_rqtl,df_pqtl)
df$qtl<-factor(df$qtl,levels=c('eQTL','rQTL','pQTL'),labels=c('eQTL','rQTL','pQTL'))
#head(df)
log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))
library(RColorBrewer)
display.brewer.all(type = "qual")
display.brewer.pal(8,"Set2")
colors = brewer.pal(8,"Set2")[1:3]
colors
p_qqplot<-ggplot(df,aes(expected, observed)) +
    geom_point(size = 3,aes(color=qtl,shape=qtl)) +
    scale_color_manual(name="QTL number",labels=c(paste0('eQTL  ',eqtl_sig_num),paste0('rQTL  ',rqtl_sig_num),paste0('pQTL  ',pqtl_sig_num)),breaks=c('eQTL','rQTL','pQTL'),values=colors)+
    scale_shape_manual(name="QTL number",labels=c(paste0('eQTL  ',eqtl_sig_num),paste0('rQTL  ',rqtl_sig_num),paste0('pQTL  ',pqtl_sig_num)),breaks=c('eQTL','rQTL','pQTL'),values=c(15,16,17))+
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +ylab(log10Po)+
    theme_bw()+
    theme(legend.position=c(0.2,0.8),legend.text=element_text(size=17),legend.title=element_text(size=15))+
    theme(axis.title.x = element_text(margin = margin(0.2,0,0,1,'cm'),size=20),axis.text.x=element_text(size=12),
        axis.title.y = element_text(margin = margin(0,0,0,0,'cm'),size=20),axis.text.y=element_text(size=12))
    #scale_fill_manual(labels=c('eQTL 2140','rQTL 840','pQTL 246'),breaks=c('eQTL','rQTL','pQTL'),values=c('green','blue','red'))+
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
############ 1b
cor<-read.table('pi1/pi1',sep="\t",header=T)
library(ggplot2)
library(ggsci)
library(ggpattern)
####
library(RColorBrewer)
display.brewer.all(type = "qual")
display.brewer.pal(8,"Set2")
colors = brewer.pal(8,"Set2")[1:3]
#library(viridis)
library(corrplot);library(ggplot2);library(ggplotify)
cor$replicate<-factor(cor$replicate,levels=c('eQTL','rQTL','pQTL'),labels=c('mRNA','ribosome\noccupancy','protein'))
cor$discovery<-factor(cor$discovery,levels=c('eQTL','rQTL','pQTL'),labels=c('mRNA','ribosome\noccupancy','protein'))
p_pi1<-ggplot(cor, aes(x=replicate,y=discovery)) +
  xlab("Replication omics types") + ylab('Discovery omics types') +
  theme_bw() + theme(panel.grid.major = element_blank()) + theme(legend.key=element_blank()) +
  theme(axis.text.x=element_text(angle=0,hjust=0.5, vjust=0.5,size=16),axis.text.y=element_text(angle=0,hjust=0.5, vjust=0.5,size=16)) +
  geom_tile(aes(fill=pi1)) + labs(fill=expression(paste(pi,"1",sep="")))+
  coord_equal()+
  geom_text(aes(x=replicate, y=discovery, label = pi1),size=8) +
  theme(axis.title.x = element_text(margin = margin(0.2,0.2,0,1,'cm'),size=20),
        axis.title.y = element_text(margin = margin(0,0.2,0,0,'cm'),size=20))+
  theme(legend.text=element_text(size=20),legend.title=element_text(size=20))+
  scale_fill_gradient2(limits=c(0,1),breaks=c(0,0.5,1),low="white",mid="#7EA772",high="#006400",midpoint=0.5)
######### 1c
library(plotrix)
library(ggplot2)
eqtl<-read.table('effectSize/eqtl.txt',sep="\t")
eqtl_mean<-mean(abs(eqtl[,3]))
eqtl_se<-std.error(abs(eqtl[,3]))
eqtl_sd<-sd(abs(eqtl[,3]))
eqtl_upper<-quantile(abs(eqtl[,3]),probs=0.25)
eqtl_lower<-quantile(abs(eqtl[,3]),probs=0.75)
eqtl_ci_upper<-CI(abs(eqtl[,3]),ci=0.95)[1]
eqtl_ci_lower<-CI(abs(eqtl[,3]),ci=0.95)[3]
#print(paste(eqtl_ci_upper,eqtl_mean,eqtl_ci_lower))
colnames(eqtl)<-c('gene','snp','slope')
eqtl$omics<-rep('eQTL',dim(eqtl)[1])

rqtl<-read.table('effectSize/rqtl.txt',sep="\t")
rqtl_mean<-mean(abs(rqtl[,3]))
rqtl_se<-std.error(abs(rqtl[,3]))
rqtl_sd<-sd(abs(rqtl[,3]))
rqtl_upper<-quantile(abs(rqtl[,3]),probs=0.25)
rqtl_lower<-quantile(abs(rqtl[,3]),probs=0.75)
rqtl_ci_upper<-CI(abs(rqtl[,3]),ci=0.95)[1]
rqtl_ci_lower<-CI(abs(rqtl[,3]),ci=0.95)[3]
colnames(rqtl)<-c('gene','snp','slope')
rqtl$omics<-rep('rQTL',dim(rqtl)[1])

pqtl<-read.table('effectSize/pqtl.txt',sep="\t")
pqtl_mean<-mean(abs(pqtl[,3]))
pqtl_se<-std.error(abs(pqtl[,3]))
pqtl_sd<-sd(abs(pqtl[,3]))
pqtl_upper<-quantile(abs(pqtl[,3]),probs=0.25)
pqtl_lower<-quantile(abs(pqtl[,3]),probs=0.75)
pqtl_ci_upper<-CI(abs(pqtl[,3]),ci=0.95)[1]
pqtl_ci_lower<-CI(abs(pqtl[,3]),ci=0.95)[3]
colnames(pqtl)<-c('gene','snp','slope')
pqtl$omics<-rep('pQTL',dim(pqtl)[1])

print(paste(eqtl_mean,rqtl_mean,pqtl_mean,sep=","))
#print(paste(eqtl_se,rqtl_se,pqtl_se,sep=","))
#print(paste(eqtl_sd,rqtl_sd,pqtl_sd,sep=","))
print(paste(eqtl_ci_upper,rqtl_ci_upper,pqtl_ci_upper,sep=","))
print(paste(eqtl_ci_lower,rqtl_ci_lower,pqtl_ci_lower,sep=","))
#wilcox.test(abs(eqtl[,3]),abs(rqtl[,3]),alternative="greater")
#wilcox.test(abs(eqtl[,3]),abs(pqtl[,3]),alternative="greater")
#wilcox.test(abs(rqtl[,3]),abs(pqtl[,3]),alternative="greater")
##
#print('two-sided')
#wilcox.test(abs(eqtl[,3]),abs(rqtl[,3]))
#wilcox.test(abs(eqtl[,3]),abs(pqtl[,3]))
#wilcox.test(abs(rqtl[,3]),abs(pqtl[,3]))
t.test(abs(eqtl[,3]),abs(rqtl[,3]),alternative="greater")
t.test(abs(rqtl[,3]),abs(pqtl[,3]),alternative="greater")
t.test(abs(eqtl[,3]),abs(pqtl[,3]),alternative="greater")


data<-rbind(eqtl,rqtl,pqtl)
data$omics<-factor(data$omics,levels=c('eQTL','rQTL','pQTL'),labels=c('eQTL','rQTL','pQTL'))
p_effect_size <- ggplot(data, aes(x=omics, y=abs(slope),color=omics)) +
    geom_boxplot()+
    labs(x='',y='Effect size')+
    scale_x_discrete(labels = c("eQTL" = "mRNA", "rQTL" = "ribosome\noccupancy","pQTL"="protein"))+
    theme_bw()+
    scale_color_manual(breaks=c('eQTL','rQTL','pQTL'),values=colors)+
    theme(legend.position = 'none')+
    theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

data2<-data.frame(
   x=factor(c('eQTL','rQTL','pQTL'),levels=c('eQTL','rQTL','pQTL'),labels=c('eQTL','rQTL','pQTL')),
   y=c(eqtl_mean,rqtl_mean,pqtl_mean),
   #sd=c(eqtl_sd,rqtl_sd,pqtl_sd))
   upper_ci=c(eqtl_ci_upper,rqtl_ci_upper,pqtl_ci_upper),
   lower_ci=c(eqtl_ci_lower,rqtl_ci_lower,pqtl_ci_lower))
p_effect_size2 <- ggplot(data2, aes(x, y,color=x)) + geom_point() +
#   geom_errorbar(aes(ymax = y + sd,
#                     ymin = y - sd))+
   geom_errorbar(aes(ymax = upper_ci,
                     ymin = lower_ci))+
   ylim(c(0,0.4))+
   labs(x='',y='Effect size')+
   scale_color_manual(breaks=c('eQTL','rQTL','pQTL'),values=colors)+
   scale_x_discrete(labels = c("eQTL" = "mRNA", "rQTL" = "ribosome\noccupancy","pQTL"="protein"))+
   theme_bw()+
   theme(legend.position='none')+
   theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
         axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

library(patchwork)
design<-'AAAAA
        AAAAA
        BB#CC
        BB#CC'
p<-wrap_plots(A=p_qqplot,B=p_pi1,C=p_effect_size2,design=design,widths=unit(c(2,7,2,2,7), c('cm', 'cm')),heights=unit(c(7,2,4,5), c('cm', 'cm')))+plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20,face='bold'))
ggsave('fig1.png',p,width=12,height=10)
