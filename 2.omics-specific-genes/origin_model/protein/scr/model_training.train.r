i=as.numeric(commandArgs(TRUE))
main_dir='/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code_demo/2.omics-specific-genes/origin_model/protein'
source(paste0(main_dir,'/scr/function.r'))

model_training(main_dir=main_dir,
               cis_window_size=1e6,
               plink_file_name<-'genotypes.sub.addCM',
               expression_file_name<-'sub.expr',
               annotation_file_name<-'gencode.v18.coding.gene.txt',
               i=i
               )

