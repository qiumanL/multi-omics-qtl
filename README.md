Code repository related to manuscript "The impact of common variants on gene expression in the human brain: from RNA to protein to schizophrenia risk."

Following the run.sh in each content to run the code.

# Requirements

## softwares

     QTLTools v1.3.1 https://qtltools.github.io/qtltools/
     GNU parallel 20131022 https://ftp.gnu.org/gnu/parallel/ (if you use job scheduling system, like PBS, slurm, just leave this software)
     MESC https://github.com/douglasyao/mesc/wiki
     PLINK v1.90p https://www.cog-genomics.org/plink2/
     GCTA v1.26.0 https://yanglab.westlake.edu.cn/software/gcta/pre_gcta/gcta_1.26.0.zip
     MR-JTI https://github.com/gamazonlab/MR-JTI/blob/master/model_training/predixcan/src/ (a part of codes used in 2.omics-specific-genes)
     MetaXcan v0.6.0 https://github.com/hakyimlab/MetaXcan (3.scz_risk_genes)
     python 3.6.3 https://www.python.org/downloads/release/python-363/
     R v4.0.3 https://cran.r-project.org/src/base/R-4/R-4.0.3.tar.gz
     libraries: 
          qvalue 2.24.0; dplyr 1.0.10; ggplot2 3.4.0; ggrepel 0.9.1; ggplotify 0.1.0; cowplot 1.1.1; patchwork 1.1.1; ggpubr 0.4.0;
          ggVennDiagram 1.2.0; RColorBrewer 1.1.2; susieR 0.12.27; ivreg 0.6.1; TwoSampleMR 0.5.6; getopt 1.20.3; 
     
     
## datasets
     gencode v19 annotation files https://www.gencodegenes.org/human/release_19.html
     GTEx BA9 eQTL https://storage.googleapis.com/adult-gtex/bulk-qtl/v7/single-tissue-cis-qtl/all_snp_gene_associations/Brain_Frontal_Cortex_BA9.allpairs.txt.gz 
     schizophrenia GWAS PGC3 https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/28169757/PGC3_SCZ_wave3_public.v2.tsv.gz?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAIYCQYOYV5JSSROOA/20240123/eu-west-1/s3/aws4_request&X-Amz-Date=20240123T213925Z&X-Amz-Expires=10&X-Amz-SignedHeaders=host&X-Amz-Signature=72df729af40cb7305ba5e23f7e0313b7378c4a21ffe09ff885779cce55c93d5b
     1000G http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/
