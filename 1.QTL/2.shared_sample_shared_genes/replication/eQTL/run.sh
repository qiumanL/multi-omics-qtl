## Shared genes and snps
time awk 'NR==FNR&&FNR>1{a[$1]++}NR>FNR&&!b[$1]++&&($1 in a){print $1}' Brain_Frontal_Cortex_BA9.allpairs.txt.format <(zcat eqtl.nominal.allpair.chr{1..22}.gz) > share.gene &
time awk 'NR==FNR&&FNR>1{a[$2]++}NR>FNR{k=$9"_"$10;if(!b[k]++&&(k in a))print k}' Brain_Frontal_Cortex_BA9.allpairs.txt.format <(zcat eqtl.nominal.allpair.chr{1..22}.gz) > share.snp &

## Extract qtls
time awk 'ARGIND==1{g[$1]++}ARGIND==2{s[$1]++}ARGIND==3&&($1 in g)&&($2 in s){print $0}' share.gene share.snp Brain_Frontal_Cortex_BA9.allpairs.txt.format |awk '!a[$1$2]++' > eqtl.gtex.nominal.chrall &
time awk 'ARGIND==1{g[$1]++}ARGIND==2{s[$1]++}ARGIND>=3{k=$9"_"$10;if(($1 in g)&&(k in s)){print $1"\t"k"\t"$12"\t"$14}}' share.gene share.snp <(zcat eqtl.nominal.allpair.chr{1..22}.gz) |awk '!a[$1$2]++' > eqtl.braingvex.nominal.chrall &

## Pi1
echo Proportion of non-null hypothesis in gtex qtl pair
for t in 0.01 0.0001 0.000001 0.00000001; do
    awk 'ARGIND==1&&$3<'$t'{a[$1$2]=1}ARGIND==2&&($1$2 in a){print $3}' eqtl.braingvex.nominal.chrall eqtl.gtex.nominal.chrall > pvalList.braingvex-gtex.$t
    Rscript qvalue.R pvalList.braingvex-gtex.$t 0.01
done
echo Proportion of non-null hypothesis in braingvex qtl pair
for t in 0.01 0.0001 0.000001 0.00000001; do
    awk 'ARGIND==1&&$3<'$t'{a[$1$2]=1}ARGIND==2&&($1$2 in a){print $3}' eqtl.gtex.nominal.chrall eqtl.braingvex.nominal.chrall > pvalList.gtex-braingvex.$t
    Rscript qvalue.R pvalList.gtex-braingvex.$t 0.01
done

echo overlap of qtl pairs
for t in 0.01 0.0001 0.000001 0.00000001; do
    awk 'ARGIND==1&&$3<'$t'{a[$1$2]=1}ARGIND==2&&$3<'$t'{if($1$2 in a) share[$1$2]=1; b[$1$2]=1}END{print length(a),length(b),length(share)}' eqtl.braingvex.nominal.chrall eqtl.gtex.nominal.chrall
done

echo overlap of qtl genes
for t in 0.01 0.0001 0.000001 0.00000001; do
    awk 'ARGIND==1&&$3<'$t'{a[$1]=1}ARGIND==2&&$3<'$t'{if($1 in a) share[$1]=1; b[$1]=1}END{print length(a),length(b),length(share)}' eqtl.braingvex.nominal.chrall eqtl.gtex.nominal.chrall
done

echo SNP-gene overlap rate: The SNP-gene overlap rates were calculated based on the percentage of shared SNPs associated with the same gene. 
for t in 0.01 0.0001 0.000001 0.00000001; do
    awk 'NR==FNR&&$3<'$t'{a[$2]=1}NR>FNR&&$3<'$t'&&($2 in a)&&!b[$2]++{print $2}' eqtl.braingvex.nominal.chrall eqtl.gtex.nominal.chrall > snp.shared.braingvex-gtex.sig
    total=$(cat snp.shared.braingvex-gtex.sig|wc -l)
    awk 'ARGIND==1{v[$1]=1}ARGIND==2&&($2 in v)&&$3<'$t'{a[$1$2]=1}ARGIND==3&&($2 in v)&&$3<'$t'{if($1$2 in a) share[$2]++}END{print length(share),'$total',length(share)/'$total'}' snp.shared.braingvex-gtex.sig eqtl.braingvex.nominal.chrall eqtl.gtex.nominal.chrall
done

