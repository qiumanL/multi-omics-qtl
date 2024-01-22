
## Revised based on https://github.com/gusevlab/fusion_twas/blob/master/FUSION.post_process.R
## Paper: Integrative approaches for large-scale transcriptome-wide association studies

## https://github.com/gusevlab/fusion_twas/issues/13
library("BEDMatrix")
read_plink_custom <- function(root, impute = c('none', 'avg', 'random')) {
    if(impute == 'random') {
        stop("The 'impute' random option has not been implemented.", call. = FALSE)
    }
    
    ## structure from https://github.com/gabraham/plink2R/blob/master/plink2R/R/plink2R.R
    proot <- path.expand(root)
    
    bedfile <- paste(proot, ".bed", sep="")
    famfile <- paste(proot, ".fam", sep="")
    bimfile <- paste(proot, ".bim", sep="")
    
    ## Could change this code to use data.table
    bim <- read.table(bimfile, header=FALSE, sep="", stringsAsFactors=FALSE)
    fam <- read.table(famfile, header=FALSE, sep="", stringsAsFactors=FALSE)
    ## Set the dimensions
    geno <- BEDMatrix::BEDMatrix(bedfile, n = nrow(fam), p = nrow(bim))
    
    ## Convert to a matrix
    geno <- as.matrix(geno)
    if(impute == 'avg') {
        ## Check if any are missing
        geno_na <- is.na(geno)
        if(any(geno_na)) {
            means <- colMeans(geno, na.rm = TRUE)
            geno[geno_na] <- rep(means, colSums(geno_na))
        }
    }
    colnames(geno) <- bim[,2]
    rownames(geno) <- paste(fam[,1], fam[, 2], sep=":")

    list(bed=geno, fam=fam, bim=bim)    
}

## read weight
cat("Reading weights...\n")
wgt<-read.table("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/protein/output/braingvex.db.weights",head=T,sep="\t",stringsAsFactors=F)
wgt<-wgt[!duplicated(paste(wgt$rsid,wgt$gene,sep=".")),]

## read genotypes
cat("Reading genotypes...\n")
genos <- read_plink_custom("/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/protein/data/genotypes.sub.addCM", impute = 'avg')
#MAFS = apply(genos$bed,2,mean)
genos$bed = scale(genos$bed)
rownames(genos$bed) <- do.call("rbind",strsplit(rownames(genos$bed),split=":"))[,1]
genos$bed<-genos$bed[,unique(wgt$rsid)]  # subset snps, to save memory and improve speed.

max_r2<-0.9  # remove splicing events with correlation larger than this value

## Start calculating
library(reshape)
#for(filename in dir(".","*.gene$")){
#for(filename in c("25056061.2014Nature.sczVScontrol.asso.csv.gene","29483656.2018NG.sczVScontrol.asso.csv.gene","29906448.2018Cell.sczVScontrol.asso.csv.gene")){
for(filename in c("2021PGC3.asso.csv.gene")){
    #filename<-"29483656.2018NG.sczVScontrol.asso.csv.gene"
    cat("gene\tnProtein\tminPval\tnProtein_noLD\tpvalue.omnibus\n", file=paste0(filename,".omnibus") , append=F )
    d<-read.table(filename,head=T,sep=",",stringsAsFactors=F)
    
    genelist<-unique(d$gene)
    nGene <- length(genelist)
    n = 1
    for(g in genelist){
        #g<-"ENSG00000234127"
        cat(paste0(filename,", #gene: ",nGene,", now: ",n,", genename: ",g,"\n"))
        n = n+1
        dd<-d[d$gene==g,]
        rownames(dd)<-dd$gene_name
        M<-nrow(dd)
        
        # read weight
        wgt.sub<-wgt[wgt$gene %in% dd$gene_name,]
        wgt.sub.cast<-cast(wgt.sub[,1:3],rsid~gene,value="weight",fill=0,drop=F)
        rownames(wgt.sub.cast)<-wgt.sub.cast$rsid
        wgt.sub.cast<-wgt.sub.cast[,-1,drop=F]
        wgt.matrix<-matrix(as.numeric(unlist(wgt.sub.cast)),nrow=nrow(wgt.sub.cast))
        rownames(wgt.matrix)<-rownames(wgt.sub.cast)
        colnames(wgt.matrix)<-colnames(wgt.sub.cast)
        
        # read genotypes
        cur.genos<-genos$bed[,rownames(wgt.matrix)]
        
        # Correlation among predicted expression of splicing events
        # cur.genos: row = samples, col = snps
        # wgt.matrix: row = snps, col = splicing events
        ge_g.matrix = cur.genos %*% wgt.matrix
        ge_g.cor = cor(ge_g.matrix)
        ge_g.cor[is.na(ge_g.cor)] <- 0
        
        if(M==1){
            cat( paste( g, M, dd$pvalue, M, dd$pvalue, sep="\t"), '\n', sep='', file=paste0(filename,".omnibus"), append=T )
            next
        }
        
        # do "informed" LD-pruning to remove highly correlated splicing events
        dd <- dd[colnames(ge_g.matrix),]
        pruned = rep(F,M)
        for(i in order(dd$zscore^2,decreasing=T) ) {
            if ( !pruned[i] ) {
                # remove anything in LD
                pruned[ ge_g.cor[ i , ]^2 > max_r2 ] = T
                pruned[i] = F
            }
        }
        
        # perform the omnibus test
        if ( sum(!pruned) > 1 ) {
            chisq = t(dd$zscore[!pruned]) %*% solve(ge_g.cor[!pruned,!pruned]) %*% dd$zscore[!pruned]
            pv.chi = pchisq( chisq , df=sum(!pruned) , lower.tail=F)
            cat( paste( g, M, min(dd$pvalue), sum(!pruned), pv.chi[1,1], sep="\t"), '\n', sep='', file=paste0(filename,".omnibus"), append=T )
        }
    }
}
