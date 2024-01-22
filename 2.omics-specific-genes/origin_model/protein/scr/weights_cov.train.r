main_dir='/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code_demo/2.omics-specific-genes/origin_model/protein'
source(paste0(main_dir,'/scr/function.r'))

weights_cov(main_dir=main_dir,
            plink_file_name<-'genotypes.sub.addCM',
            expression_file_name<-'sub.expr',
            output_file_name<-'braingvex'
)



