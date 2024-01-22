main_dir='/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/conditional_model/rna'

source(paste0(main_dir,'/scr/function.r'))

weights_cov(main_dir=main_dir,
            plink_file_name<-'genotypes.sub.addCM',
            expression_file_name<-'sub.expr.resid',
            output_file_name<-'braingvex'
)



