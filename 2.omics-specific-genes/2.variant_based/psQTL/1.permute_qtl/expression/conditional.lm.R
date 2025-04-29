options(stringsAsFactors=FALSE)
argv<-commandArgs()
## Initialize
typelist_all<-c("eQTL","rQTL","pQTL")
typeObs<-as.character(argv[6])
typelist<-typelist_all[-which(typelist_all==typeObs)]
typelist

## read observed and predicted values
# for protein sub expr, ensure the first column contains only gene ids
observed<-as.matrix(read.table(paste0("1.QTL/shared_sample_shared_genes/",typeObs,"/expression/sub.expr"),head=T,row.names=1,sep="\t",check.names=F))
predicted<-list()
for(type in typelist){
    predicted[[type]]<-as.matrix(read.table(paste0("1.QTL/shared_sample_shared_genes/",type,"/expression/sub.expr"),head=T,row.names=1,sep="\t",check.names=F))
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
# rna / ribo / protein
genelist<-rownames(read.table("1.QTL/shared_sample_shared_genes/eQTL/expression/sub.expr",head=T,row.names=1,sep="\t"))
iso2gene[["rna"]]<-genelist
iso2gene[["ribo"]]<-genelist
iso2gene[["protein"]]<-genelist
names(iso2gene[["rna"]])<-genelist
names(iso2gene[["ribo"]])<-genelist
names(iso2gene[["protein"]])<-genelist
for(g in genelist){
    gene2iso[["rna"]][[g]]<-g
    gene2iso[["ribo"]][[g]]<-g
    gene2iso[["protein"]][[g]]<-g
}

## re-construct predicted values by geneList
have_removed_gene<-vector()
covs<-list()
for(g in iso2gene[[typeObs]][rownames(observed)]){
    d<-matrix(1,nrow=ncol(observed),ncol=1)
    colnames(d)<-"intercept"
    for(type in typelist){
            for(i0 in gene2iso[[type]][[g]]){
                if(i0 %in% rownames(predicted[[type]])){
                    d<-cbind(d,t(predicted[[type]][i0,,drop=F]))
                    colnames(d)<-c(colnames(d)[1:(ncol(d)-1)],paste(type,i0,sep="."))
            }
        }
    }
    covs[[g]]<-d
}
##
g<-"ENSG00000004700"
head(covs[[g]])
nrow(observed)

## get residuals
mod.sum<-summary(lm(observed[g,]~covs[[iso2gene[[typeObs]][g]]]))
names(mod.sum)
mod.sum$r.squared
r2_data<-as.data.frame(sapply(rownames(observed),function(g){summary(lm(observed[g,]~covs[[iso2gene[[typeObs]][g]]]))$r.squared}))
colnames(r2_data)<-'R2'
head(r2_data)
write.table(r2_data,"remove.r2.txt", sep="\t", row.names=T, quote=F,col.names=T)

resids<-t(sapply(rownames(observed),function(g){residuals(lm(observed[g,]~covs[[iso2gene[[typeObs]][g]]]))+mean(observed[g,])}))
write.table(resids, file="sub.expr.resid", sep="\t", row.names=T, quote=F,col.names=NA)

#save.image("conditional.lm.RData")
