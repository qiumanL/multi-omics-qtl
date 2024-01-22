
model_training<-function(main_dir,
                         plink_file_name,
                         expression_file_name,
                         annotation_file_name,
                         cis_window_size=1e6,
                         i){
  library(glmnet)
  # main_dir<-'/gpfs23/data/coxvgi/zhoud2/projects/predixcan'
  # cis_window_size=1e6
  # plink_file_name<-'geuvadis421'
  # expression_file_name<-'geuvadis_residual.txt'
  # annotation_file_name<-'gencode.v18.coding.gene.txt'
  # i=1
  
  #mkdir
  cmd<-paste0('mkdir ',main_dir,'/gene/; mkdir ',main_dir,'/weights/; mkdir ',main_dir,'/cov/; mkdir ',main_dir,'/output/')
  system(cmd,wait = T)
  
  #load annotation
  anno<-read.table(paste0(main_dir,'/data/',annotation_file_name),header = T,stringsAsFactors = F) #keeps the colnames consistant with the exsample file
  #load expression data
  print('loading expression data...')
  exp_all.h<-as.character(read.table(paste0(main_dir,'/data/',expression_file_name),nrow=1,header = F,stringsAsFactors = F))
  exp_all<-read.table(paste0(main_dir,'/data/',expression_file_name),skip=(i-1)*100+1,nrow=100,header = F,stringsAsFactors = F)
  colnames(exp_all)<-exp_all.h
  #exp_all<-read.table(paste0(main_dir,'/data/',expression_file_name),header = T,stringsAsFactors = F)  #keeps IID in genotype file match with IID in expression data
  #gene_list
  gene_list<-intersect(exp_all$geneid,anno$geneid)
  
  #i_start<-(i-1)*100+1
  #i_end<-min(i*100,length(gene_list))
  i_start<-1
  i_end<-length(gene_list)
  
  if (i_start<=length(gene_list)){
    print(paste0('from No.',i_start,' gene to No.',i_end,' gene...'))
    
    for (i in i_start:i_end){
      geneid=gene_list[i] #geneid='ENSG00000164978'
      print(paste0('processing No.',i,' ',geneid))
      
      #get gene position
      chr=sub('^...','',anno$chr[which(anno$geneid==geneid)])
      pos_l=max(anno$pos_l[which(anno$geneid==geneid)]-cis_window_size,1)
      pos_r=anno$pos_r[which(anno$geneid==geneid)]+cis_window_size
      
      #generate expression df
      exp<-data.frame('IID'=as.character(colnames(exp_all)[-1]),'exp'=as.numeric(exp_all[exp_pos<-which(exp_all[,1]==geneid),][-1]))
      
      each_gene(main_dir=main_dir,plink_file_name=plink_file_name,chr=chr,pos_l=pos_l,pos_r=pos_r,geneid=geneid,exp=exp,iGene=i)
      
    }
  }else{
    print('start_i > N of genes')
  }

}

#----------------------------------------


each_gene<-function(main_dir,plink_file_name,chr,pos_l,pos_r,geneid,exp,iGene){
  #generate genotype dosage file for each gene
  #extract from plink
  cmd<-paste0('plink --bfile ',main_dir,'/data/',plink_file_name,' --chr ',chr,' --from-bp ',pos_l,' --to-bp ',pos_r,' --recodeA --out ',main_dir,'/gene/',geneid)
  system(cmd,wait = T)
  
  #load dosage
  dosage_raw<-try(read.table(paste0(main_dir,'/gene/',geneid,'.raw'),header = T,stringsAsFactors = F))
  if(!('try-error' %in% class(dosage_raw))){
    #rm raw file
    cmd<-paste0('rm ',main_dir,'/gene/',geneid,'*'); system(cmd,wait = T)
    #del useless cols
    dosage<-dosage_raw[,-c(1,3:6)] 
    
    #merge genotype and expression
    df<-merge(exp,dosage,by='IID')
    
    #run elastic net
    y<-as.matrix(df[,2])
    x<-as.matrix(df[,c(3:ncol(df))])
    rm(df)
    #set.seed(as.numeric(sub('^....','',geneid)))
    set.seed(iGene)
    fit<-cv.glmnet(x=x,y=y, nfolds = 10,keep = T,alpha=0.5)
    beta=as.numeric(fit$glmnet.fit$beta[,which(fit$lambda==fit$lambda.min)])
    
    performance<-cor.test(y,fit$fit.preval[,which(fit$lambda == fit$lambda.min)])
    r<-as.numeric(performance$estimate)
    p<-as.numeric(performance$p.value)
    
    #generate weight df
    weights_df<-data.frame('rsid'=colnames(dosage)[-1],'weights'=beta,stringsAsFactors = F,'r2'=r^2,'p'=p)
    weights_df$effect_allele<-sapply(weights_df$rsid, function(x) strsplit(x,"[_]")[[1]][2])
    weights_df$rsid<-sapply(weights_df$rsid, function(x) strsplit(x,"[_]")[[1]][1])
    #rm weight=0 rows
    weights_df<-weights_df[weights_df$weights!=0,]
    
    if (nrow(weights_df)>0 & r>0.1 & p<0.05){
      
      #generate covariance matrix
      colnames(dosage)[-1]<-sapply(colnames(dosage)[-1], function(x) strsplit(x,"[_]")[[1]][1])
      dosage<-dosage[,c(which(colnames(dosage) %in% weights_df$rsid)),drop=FALSE]
      cov_df<-as.data.frame(matrix(data=NA,ncol=4,nrow=0)) 
      colnames(cov_df)<-c('GENE','RSID1','RSID2','VALUE')
      #calculate covariance
      snp_list<-colnames(dosage)
      o_i=1
      for (k in 1:length(snp_list)){
        for (j in k:length(snp_list)){
          cov_df[o_i,2]<-snp_list[k]
          cov_df[o_i,3]<-snp_list[j]
          cov_df[o_i,4]<-round(cov(dosage[,which(colnames(dosage)==snp_list[k])],dosage[,which(colnames(dosage)==snp_list[j])]),3)
          o_i=o_i+1
        }
      }
      cov_df[,1]<-geneid
      
      #output cov
      write.table(cov_df,paste0(main_dir,'/cov/',geneid,'.txt'),sep='\t',quote = F,row.names = F)
      #output weights
      write.table(weights_df,paste0(main_dir,'/weights/',geneid,'.txt'),sep='\t',quote = F,row.names = F)
    }
  }else{
    print('no genotype data for this gene')
  }
}

#------------------------------------------------

weights_cov<-function(main_dir,plink_file_name,expression_file_name,output_file_name){
  
  #gene list
  gene_list<-sub('....$','',dir(paste0(main_dir,'/cov/')))
  
  #load weights
  print('loading weights...')
  db_weights<-list()
  n_rows=0
  for (i in 1:length(gene_list)){
    geneid<-gene_list[i]
    db_weights[[i]]<-read.table(paste0(main_dir,'/weights/',geneid,'.txt'),header = T,stringsAsFactors = F,colClasses=c("character",rep("numeric",3),"character"))
    n_rows=n_rows+nrow(db_weights[[i]])
  }
  
  #db.weight
  weights<-as.data.frame(matrix(data=NA,ncol=4,nrow=n_rows))
  colnames(weights)<-c('rsid','gene','weight','eff_allele')
  
  #db.construction
  construction<-as.data.frame(matrix(data=NA,ncol=2,nrow=length(gene_list)))
  colnames(construction)<-c('chr','cv.seed')
  construction$chr=1 #it doesn't matter
  
  #db.sample_info
  n_sample<-read.table(paste0(main_dir,'/data/',expression_file_name),header = T,nrows=1,stringsAsFactors = F)
  sample_info<-data.frame('n.samples'=ncol(n_sample)-1)
  
  #db.extra
  extra<-as.data.frame(matrix(data=NA,ncol=6,nrow=length(gene_list)))
  colnames(extra)<-c('gene','genename','pred.perf.R2','n.snps.in.model','pred.perf.pval','pred.perf.qval')
  
  n_rows=0
  for (i in 1:length(gene_list)){
    #construction.seed
    construction[i,2]<-as.numeric(sub('^....','',gene_list[i]))
    #extra.gene
    extra[i,1:2]<-gene_list[i]
    #extra.r2
    extra[i,3]<-db_weights[[i]][1,3]
    #extra.n.snps
    extra[i,4]<-nrow(db_weights[[i]])
    #extra.n.snps
    extra[i,5:6]<-db_weights[[i]][1,4] #the qval is not valid
    
    weights[(n_rows+1):(n_rows+nrow(db_weights[[i]])),]<-as.matrix(data.frame('rsid'=db_weights[[i]]$rsid,'gene'=gene_list[i],'weight'=db_weights[[i]]$weights,'effect_allele'=db_weights[[i]]$effect_allele))
    n_rows=n_rows+nrow(db_weights[[i]])
  }
  
  #load bim to get ref allele
  bim<-read.table(paste0(main_dir,'/data/',plink_file_name,'.bim'),header = F,stringsAsFactors = F)
  bim<-bim[bim[,2]!='.',];bim<-bim[!duplicated(bim[,2]),]
  bim<-bim[which(bim[,2] %in% weights[,1]),-1]
  
  weights<-merge(weights,bim,by=1)
  weights$ref_allele<-ifelse(weights$eff_allele==weights$V5,weights$V6,weights$V5)
  weights<-weights[,c(1,2,3,9,4)]
  
  #output db file
  print('generating weights into .db...')
  library(RSQLite)
  com<-paste0("rm ",main_dir,"/output/",output_file_name,".db");system(com,wait = T)
  
  db<-dbConnect(RSQLite::SQLite(), paste0(main_dir,"/output/",output_file_name,".db"))
  dbWriteTable(db, "weights", weights)
  dbWriteTable(db, "construction", construction)
  dbWriteTable(db, "extra", extra)
  dbWriteTable(db, "sample_info", sample_info)
  dbDisconnect(db)
  
  #load cov
  print('loading cov...')
  cov<-list()
  n_rows=0
  for (i in 1:length(gene_list)){
    geneid<-gene_list[i]
    cov[[i]]<-read.table(paste0(main_dir,'/cov/',geneid,'.txt'),header = T,stringsAsFactors = F)
    n_rows=n_rows+nrow(cov[[i]])
  }
  cov_df<-as.data.frame(matrix(data=NA,ncol=4,nrow=n_rows)) 
  colnames(cov_df)<-c('GENE','RSID1','RSID2','VALUE')
  #combine cov
  print('combine cov...')
  n_rows=0
  for (i in 1:length(gene_list)){
    cov_df[(n_rows+1):(n_rows+nrow(cov[[i]])),]<-cov[[i]]
    n_rows=n_rows+nrow(cov[[i]])
  }
  
  #output cov
  write.table(cov_df,paste0(main_dir,'/output/',output_file_name,'.cov'),sep='\t',row.names = F,quote = F)
  print('finished...')
}


