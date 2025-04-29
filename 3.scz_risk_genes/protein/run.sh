python3 MetaXcan-0.6.11/software/MetaXcan.py --model_db_path braingvex.db --covariance braingvex.cov --gwas_file 2021_PGC3_SCZ.gz --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --or_column OR --pvalue_column P --output_file 2021PGC3.asso.csv

ls *.csv |while read i; do awk 'BEGIN{FS="[\t,]";OFS=","}NR==FNR{a[$1]=$3}NR>FNR{if(FNR==1){print}else if($1 in a){$1=a[$1];print $0}}' gene2prot $i > $i.gene; done

Rscript omnibus.R

for i in *.csv.gene.omnibus; do Rscript multipleTesting.omnibus.R $i; done

ls *.omnibus.multiTest |while read i; do cat $i |awk '$7<0.05{print $1}' > $i.bonferroni005; done

for i in *.bonferroni005; do awk 'NR==FNR{split($4,s,".");a[s[1]]=1}NR>FNR&&!($1 in a)' gene.MHC.bed $i > $i.nonMHC; done

awk -F "\t" 'ARGIND==1{split($1,x,"chr");a[$7]=x[2];b[$7]=($2+0+$3+0)/2;c[$7]=$6}ARGIND==2{if($1 in a){print $1"\t"a[$1]"\t"b[$1]"\t"$7"\t"c[$1]}}' gene2coor 2021PGC3.asso.csv.gene.omnibus.multiTest >2021PGC3.asso.csv.gene.omnibus.multiTest.mahattan

awk -F "\t" 'ARGIND==1{a[$7]=$6}ARGIND==2{if($1 in a){print $1"\t"a[$1]}}' gene2coor 2021PGC3.asso.csv.coding.multiTest.bonferroni005.nonMHC >2021PGC3.asso.csv.coding.multiTest.bonferroni005.nonMHC.genename

