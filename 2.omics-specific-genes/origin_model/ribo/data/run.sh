 cat /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/1.QTL/shared_sample_shared_genes/rQTL/expression/sub.expr |sed -e'1s/^gene/geneid/' -e'1s/\t/\tX/g' -e'1s/-/./g' > sub.expr

 cat /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/1.QTL/shared_sample_shared_genes/genotype/genotypes.sub.addCM.fam |sed -e's/^/X/' -e's/ / X/' -e's/-/./' -e's/-/./' > genotypes.sub.addCM.fam

##
awk -F "\t" 'NR==1{for(i=2;i<=NF;i++){print $i}}' sub.expr >sample.order
parallel -j 6 bash subVcf.dosage.sh {} ::: {1..22}
awk -F "\t" '{print $1"\t"$1}' sample.order >sample.order.fmt

for i in {1..22};do awk 'ARGIND==1{a[$2]=1}ARGIND==2{if($1":"$2 in a) b[$1":"$2]=$6}ARGIND==3{if($2 in b){$2=b[$2]}print $0}' genotypes.sub.addCM.chr"$i".dosage /zs32/data-analysis/reflib/annovar_humandb/hg19_avsnp150.txt genotypes.sub.addCM.chr"$i".dosage |gzip > genotypes.sub.addCM.chr"$i".dosage.txt.gz;done

