## chr=3 gene=ENSG00000114904 name=NEK4
chr=$1
gene=$2
name=$3
mkdir $gene
cd $gene
################ for manhattan plot
start=`grep "$gene" /zs32/data-analysis/liucy_group/jiangyi/database/gencode.v19/gene2coor|awk -F "\t" '{print $2}'`
end=`grep "$gene" /zs32/data-analysis/liucy_group/jiangyi/database/gencode.v19/gene2coor|awk -F "\t" '{print $3}'`
awk -F "\t" '{split($1,a,"_");if(a[2]+0>='"$start"'+0-1000000 && a[2]+0<='"$end"'+0+1000000){print a[1]"\t"a[2]"\t"$1"\t"a[3]"\t"a[4]"\t"$7"\t"$5"\t"$6"\t."}}' /zs32/data-analysis/liucy_group/liangqiuman/psychENCODE/shared_multi-omics/nominal_qtl/gwas/2021PGC3/2021PGC3.sczVScontrol.ma.chr"$chr" >gwas."$gene".txt
sed -i '1i#chrom\tpos\tsnp\tref\talt\tpvalue\tbeta\tstderr_beta\talt_allele_freq' gwas."$gene".txt

awk -F "\t" '{split($1,a,"_");print a[1]"\t"a[2]"\t"$1"\t"a[3]"\t"a[4]"\t"$2"\t"$3"\t"$4"\t."}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/eqtl/for_mr/chr"$chr"/eqtl.allP.forgwas.chr"$chr"."$gene" >eqtl."$gene".txt
sed -i '1i#chrom\tpos\tsnp\tref\talt\tpvalue\tbeta\tstderr_beta\talt_allele_freq' eqtl."$gene".txt

awk -F "\t" '{split($2,a,"_");print a[1]"\t"a[2]"\t"$1"\t"a[3]"\t"a[4]"\t"$3"\t"$4"\t"$5"\t."}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/rqtl/for_mr/chr"$chr"/rqtl.allP.forgwas.chr"$chr"."$gene" >rqtl."$gene".txt
sed -i '1i#chrom\tpos\tsnp\tref\talt\tpvalue\tbeta\tstderr_beta\talt_allele_freq' rqtl."$gene".txt

awk -F "\t" '{split($2,a,"_");print a[1]"\t"a[2]"\t"$2"\t"a[3]"\t"a[4]"\t"$3"\t"$4"\t"$5"\t."}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/pqtl/for_mr/chr"$chr"/pqtl.allP.forgwas.chr"$chr"."$gene" >pqtl."$gene".txt
sed -i '1i#chrom\tpos\tsnp\tref\talt\tpvalue\tbeta\tstderr_beta\talt_allele_freq' pqtl."$gene".txt
### coloc snp
top_coloc_snp=`awk -F "\t" 'NR>1{print $1}' <(head -2 /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/2.predixcan/coloc/eqtl/chr"$chr"/"$gene"/coloc_snp_PP4.txt)`
echo $top_coloc_snp >top_coloc_snp

############## for boxplot
top_eqtl_snp=`awk '$1=="'$gene'"{print $8}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/eqtl/eqtl.permute.qvalue.sig.txt`
echo $top_eqtl_snp >top_eqtl_snp
zcat /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/genotype/genotypes.all.chr"$chr".realign.vcf.gz|awk '$0~/^#CHR/{print}'|awk -F "\t" '{for(i=9;i<=NF;i++){printf $i"\n"}}' >genotype.head
zgrep "$top_eqtl_snp" /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/genotype/genotypes.all.chr"$chr".realign.vcf.gz|awk -F "\t" '{for(i=9;i<=NF;i++){printf $i"\n"}}' >genotype
paste -d "\t" genotype.head genotype >genotype.txt
### transfer to dosage
awk -F "\t" 'BEGIN{OFS="\t"}{if($2=="0|0"){$2=0;print}if($2=="0|1"||$2=="1|0"){$2=1;print}if($2=="1|1"){$2=2;print}}' genotype.txt >genotype.dosage.txt

##### RNA-Seq
###
zcat /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/eqtl/expression/sub.expr.chr"$chr".bed.gz|head -1|awk -F "\t" '{for(i=7;i<=NF;i++){printf $i"\n"}}' >exp.head.eqtl
zgrep "$gene" /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/eqtl/expression/sub.expr.chr"$chr".bed.gz|awk -F "\t" '{for(i=7;i<=NF;i++){printf $i"\n"}}' >exp.eqtl
paste -d "\t" exp.head.eqtl exp.eqtl >exp.eqtl.txt
### merge
awk -F "\t" 'ARGIND==1{a[$1]=$2}ARGIND==2{if($1 in a){print $0"\t"a[$1]"\tmRNA"}}' exp.eqtl.txt genotype.dosage.txt >all.eqtl.txt

##### Ribo-Seq
zcat /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/rqtl/expression/sub.expr.chr"$chr".bed.gz|head -1|awk -F "\t" '{for(i=7;i<=NF;i++){printf $i"\n"}}' >exp.head.rqtl
zgrep "$gene" /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/rqtl/expression/sub.expr.chr"$chr".bed.gz|awk -F "\t" '{for(i=7;i<=NF;i++){printf $i"\n"}}' >exp.rqtl
paste -d "\t" exp.head.rqtl exp.rqtl >exp.rqtl.txt
### merge
awk -F "\t" 'ARGIND==1{a[$1]=$2}ARGIND==2{if($1 in a){print $0"\t"a[$1]"\tRibo"}}' exp.rqtl.txt genotype.dosage.txt >all.rqtl.txt

##### Mass-Spec
zcat /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/pqtl/expression/sub.expr.chr"$chr".bed.gz|head -1|awk -F "\t" '{for(i=7;i<=NF;i++){printf $i"\n"}}' >exp.head.pqtl
zgrep "$gene" /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/pqtl/expression/sub.expr.chr"$chr".bed.gz|awk -F "\t" '{for(i=7;i<=NF;i++){printf $i"\n"}}' >exp.pqtl
paste -d "\t" exp.head.pqtl exp.pqtl >exp.pqtl.txt
### merge
awk -F "\t" 'ARGIND==1{a[$1]=$2}ARGIND==2{if($1 in a){print $0"\t"a[$1]"\tProtein"}}' exp.pqtl.txt genotype.dosage.txt >all.pqtl.txt

#####
cat all.eqtl.txt all.rqtl.txt all.pqtl.txt >all.txt
sed -i '1isample\tgenotype\texpression\tqtl' all.txt
