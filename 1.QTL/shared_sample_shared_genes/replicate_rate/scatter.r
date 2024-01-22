library(ggplot2)
up<-read.table('upstream/result',sep="\t",head=F)
down<-read.table('downstream/result',sep="\t",head=F)
colnames(up)<-c('discovery','replicate','cutoff','sig_overlap','overlap','rate')
colnames(down)<-c('discovery','replicate','cutoff','sig_overlap','overlap','rate')
data_merge<-rbind(up,down)
data<-subset(data_merge,discovery!="rqtl" & cutoff!=0.03 & cutoff!=0.02 & cutoff!=0.01 & cutoff!=0.008 & cutoff!=0.005 & cutoff!=0.001)
data
data$replicate<-factor(data$replicate,levels=c('eqtl','rqtl','pqtl'),labels=c("mRNA","ribosome occupancy","protein"))
data$discovery<-factor(data$discovery,levels=c('eqtl','pqtl'),labels=c("eQTL","pQTL"))
#png('replicate.png',width)
library(RColorBrewer)
display.brewer.all(type = "qual")
display.brewer.pal(8,"Set2")
colors = brewer.pal(8,"Set2")[1:3]
colors
p<-ggplot(data=data, aes(x=cutoff, y=rate,color=replicate)) + geom_point(aes(shape=replicate),size=3)+
  facet_grid(.~discovery)+
  ylab('Replication rate (%)')+
  scale_color_manual(values = colors)+
  guides(colour = guide_legend("Replication datatype"), 
           shape = guide_legend("Replication datatype"))+
  theme_bw()+ylim(0,100)+
  theme(strip.text.x = element_text(size = 20))+
  labs(x='FDR cutoff',color='Replication datatype')+
  theme(legend.position='bottom',legend.title=element_text(size=16),legend.text=element_text(size=16))+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))
ggsave('replicate.png',p,width=10,height=5)
#dev.off()
