library(ggplot2)
library(ggpubr)
data<-read.table('top10.qvalue.txt',head=T,sep="\t")
p3<-ggplot(data,aes(x=slope.185sample,y=slope.10sample))+
    geom_point()+
    theme_bw()+
    theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
         axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+
    labs(x='Effect size from 185 samples',y='Effect size from 10 samples')+
    stat_cor(method="spearman")
ggsave('cor.3copy.2fc.png',p3,height=5,width=5)

