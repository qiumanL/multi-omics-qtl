python3 /zs32/data-analysis/liucy_group/liangqiuman/software/MetaXcan/MetaXcan-0.6.11/software/MetaXcan.py --model_db_path /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/rna/output/braingvex.db --covariance /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/rna/output/braingvex.cov --gwas_file /vg_sklmgdata_hw_01/data/liangqiuman/database/PGC3/2021_PGC3_SCZ.gz --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --or_column OR --pvalue_column P --output_file 2021PGC3.asso.csv

for i in *.csv; do awk 'BEGIN{FS="[,\t]"}NR==FNR&&$5=="protein_coding"{split($4,s,".");a[s[1]]=1}NR>FNR&&(FNR==1||($1 in a)){print $0}' /zs32/data-analysis/liucy_group/jiangyi/database/gencode.v19/gencode.v19.annotation.gene.bed $i > $i.coding; done

for i in *.csv.coding; do Rscript multipleTesting.R $i; done

ls *.multiTest|while read i; do cat $i |awk '$15<0.05{print $1}' > $i.bonferroni005; done

for i in *.bonferroni005; do awk 'NR==FNR{split($4,s,".");a[s[1]]=1}NR>FNR&&!($1 in a)' /zs32/data-analysis/liucy_group/jiangyi/database/gencode.v19/gene.MHC.bed $i > $i.nonMHC; done

awk -F "\t" 'ARGIND==1{split($1,x,"chr");a[$7]=x[2];b[$7]=($2+0+$3+0)/2;c[$7]=$6}ARGIND==2{if($1 in a){print $1"\t"a[$1]"\t"b[$1]"\t"$15"\t"c[$1]}}' /zs32/data-analysis/liucy_group/jiangyi/database/gencode.v19/gene2coor 2021PGC3.asso.csv.coding.multiTest >2021PGC3.asso.csv.coding.multiTest.mahattan

awk -F "\t" 'ARGIND==1{a[$7]=$6}ARGIND==2{if($1 in a){print $1"\t"a[$1]}}' /zs32/data-analysis/liucy_group/jiangyi/database/gencode.v19/gene2coor 2021PGC3.asso.csv.coding.multiTest.bonferroni005.nonMHC >2021PGC3.asso.csv.coding.multiTest.bonferroni005.nonMHC.genename
