options(stringsAsFactors=FALSE)
argv<-commandArgs()
## Initialize
typelist_all<-c("rna","ribo","protein")
typeObs<-as.character(argv[6])
typelist<-typelist_all[-which(typelist_all==typeObs)]
typelist

## read observed and predicted values
#observed<-as.matrix(read.table(paste0("/zs32/data-analysis/liucy_group/jiangyi/psychENCODE/predixcan/subSample.rna-ribo-protein-splicing.norm2/divide2/",typeObs,"/data/expr.sub.test"),head=T,row.names=1,sep="\t"))
observed<-as.matrix(read.table(paste0("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/",typeObs,"/data/sub.expr"),head=T,row.names=1,sep="\t"))
predicted<-list()
for(type in typelist){
#    predicted[[type]]<-t(read.table(paste0("/zs32/data-analysis/liucy_group/jiangyi/psychENCODE/predixcan/subSample.rna-ribo-protein-splicing.norm2/divide2/",type,"/output/braingvex.test.predict_predicted_expression.txt2"),head=T,row.names=1,sep="\t"))
    predicted[[type]]<-t(read.table(paste0("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/",type,"/output/braingvex.predict_predicted_expression.txt2"),head=T,row.names=1,sep="\t",check.names = FALSE))
    if(any(colnames(predicted[[type]])!=colnames(observed))){
        print(paste0("colnames not matching between observed and ",type," predicted!"))
        exit(1)
    }
}

## read gene translation tables
iso2gene<-list()
gene2iso<-list()
for(type in c(typelist,typeObs)){
    gene2iso[[type]]<-list()
}
# rna / ribo
genelist<-rownames(read.table("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/rna/data/sub.expr",head=T,row.names=1,sep="\t"))
iso2gene[["rna"]]<-genelist
iso2gene[["ribo"]]<-genelist
names(iso2gene[["rna"]])<-genelist
names(iso2gene[["ribo"]])<-genelist
for(g in genelist){
    gene2iso[["rna"]][[g]]<-g
    gene2iso[["ribo"]][[g]]<-g
}
# protein
d<-read.table("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/protein/data/protein.pos.ensg",head=T,sep="\t")
iso2gene[["protein"]]<-d$ensg
names(iso2gene[["protein"]])<-d$geneid
for(g in unique(d$ensg)){
    gene2iso[["protein"]][[g]]<-d[d$ensg==g,"geneid"]
}
## re-construct predicted values by geneList
have_removed_gene<-vector()
covs<-list()
for(g in iso2gene[[typeObs]][rownames(observed)]){
    d<-matrix(1,nrow=ncol(observed),ncol=1)
    colnames(d)<-"intercept"
    for(type in typelist){
        if(g %in% iso2gene[[type]][rownames(predicted[[type]])]){
            for(i0 in gene2iso[[type]][[g]]){
                if(i0 %in% rownames(predicted[[type]])){
                    d<-cbind(d,t(predicted[[type]][i0,,drop=F]))
                    colnames(d)<-c(colnames(d)[1:(ncol(d)-1)],paste(type,i0,sep="."))
                    have_removed_gene<-append(have_removed_gene,g)
                }
            }
        }
    }
    covs[[g]]<-d
}
##
g<-"ENSG00000004700"
head(covs[[g]])
observed<-observed[iso2gene[[typeObs]][rownames(observed)] %in% have_removed_gene,]
nrow(observed)
## get residuals
resids<-t(sapply(rownames(observed),function(g){residuals(lm(observed[g,]~covs[[iso2gene[[typeObs]][g]]]))+mean(observed[g,])}))
write.table(resids, file="sub.expr.resid", sep="\t", row.names=T, quote=F,col.names=NA)

#save.image("conditional.lm.RData")
