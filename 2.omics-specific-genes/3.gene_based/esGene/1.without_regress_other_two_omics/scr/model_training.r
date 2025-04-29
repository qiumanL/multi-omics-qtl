i=as.numeric(commandArgs(TRUE))
### requiring main_dir include genotypes.sub.addCM (plink formate),sub.expr.resid and gencode.v18.coding.gene.txt
main_dir='/1.without_regress_other_two_omics'

source(paste0(main_dir,'/scr/function.r'))

model_training(main_dir=main_dir,
               cis_window_size=1e6,
               plink_file_name<-'genotypes.sub.addCM',
               expression_file_name<-'sub.expr.resid',
               annotation_file_name<-'gencode.v18.coding.gene.txt',
               i=i
               )

