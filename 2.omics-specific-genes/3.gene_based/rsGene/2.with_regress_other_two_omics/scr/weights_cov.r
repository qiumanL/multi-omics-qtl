main_dir='/2.with_regress_other_two_omics'

source(paste0(main_dir,'/scr/function.r'))

weights_cov(main_dir=main_dir,
            plink_file_name<-'genotypes.sub.addCM',
            expression_file_name<-'sub.expr.resid',
            output_file_name<-'braingvex'
)



