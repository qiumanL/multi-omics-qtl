## Shared genes and snps
time awk 'NR==FNR&&FNR>1{a[$5]++}NR>FNR{split($1,x,"|");if(!b[x[2]]++ && (x[2] in a)){print x[2]}}' rosmap.pqtl.all <(zcat pqtl.nominal.allpair.chr{1..22}.gz) > share.gene &
time awk 'NR==FNR&&FNR>1{a[$1"_"$2"_"$3"_"$4]++}NR>FNR{k=$8;if(!b[k]++&&(k in a))print k}' rosmap.pqtl.all <(zcat pqtl.nominal.allpair.chr{1..22}.gz) > share.snp &

## Extract qtls
time awk 'ARGIND==1{g[$1]++}ARGIND==2{s[$1]++}ARGIND==3{if($5 in g && $1"_"$2"_"$3"_"$4 in s){print $5"\t"$1"_"$2"_"$3"_"$4"\t"$8"\t"$6}}' share.gene share.snp rosmap.pqtl.all |awk '!a[$1$2]++' > pqtl.rosmap.nominal.chrall &
time awk 'ARGIND==1{g[$1]++}ARGIND==2{s[$1]++}ARGIND>=3{k=$8;split($1,x,"|");if((x[2] in g)&&(k in s)){print x[2]"\t"k"\t"$12"\t"$14}}' share.gene share.snp <(zcat pqtl.nominal.allpair.chr{1..22}.gz) |awk '!a[$1$2]++' > pqtl.braingvex.nominal.chrall &

## Pi1
echo Proportion of non-null hypothesis in rosmap qtl pair
for t in 0.01 0.0001 0.000001 0.00000001; do
    awk 'ARGIND==1&&$3<'$t'{a[$1$2]=1}ARGIND==2&&($1$2 in a){print $3}' pqtl.braingvex.nominal.chrall pqtl.rosmap.nominal.chrall > pvalList.braingvex-rosmap.$t
    Rscript qvalue.R pvalList.braingvex-rosmap.$t 0.01
done
echo Proportion of non-null hypothesis in braingvex qtl pair
for t in 0.01 0.0001 0.000001 0.00000001; do
    awk 'ARGIND==1&&$3<'$t'{a[$1$2]=1}ARGIND==2&&($1$2 in a){print $3}' pqtl.rosmap.nominal.chrall pqtl.braingvex.nominal.chrall > pvalList.rosmap-braingvex.$t
    Rscript qvalue.R pvalList.rosmap-braingvex.$t 0.01
done

echo overlap of qtl pairs
for t in 0.01 0.0001 0.000001 0.00000001; do
    awk 'ARGIND==1&&$3<'$t'{a[$1$2]=1}ARGIND==2&&$3<'$t'{if($1$2 in a) share[$1$2]=1; b[$1$2]=1}END{print length(a),length(b),length(share)}' pqtl.braingvex.nominal.chrall pqtl.rosmap.nominal.chrall
done

echo overlap of qtl genes
for t in 0.01 0.0001 0.000001 0.00000001; do
    awk 'ARGIND==1&&$3<'$t'{a[$1]=1}ARGIND==2&&$3<'$t'{if($1 in a) share[$1]=1; b[$1]=1}END{print length(a),length(b),length(share)}' pqtl.braingvex.nominal.chrall pqtl.rosmap.nominal.chrall
done

echo SNP-gene overlap rate: The SNP-gene overlap rates were calculated based on the percentage of shared SNPs associated with the same gene. 
for t in 0.01 0.0001 0.000001 0.00000001; do
    awk 'NR==FNR&&$3<'$t'{a[$2]=1}NR>FNR&&$3<'$t'&&($2 in a)&&!b[$2]++{print $2}' pqtl.braingvex.nominal.chrall pqtl.rosmap.nominal.chrall > snp.shared.braingvex-rosmap.sig
    total=$(cat snp.shared.braingvex-rosmap.sig|wc -l)
    awk 'ARGIND==1{v[$1]=1}ARGIND==2&&($2 in v)&&$3<'$t'{a[$1$2]=1}ARGIND==3&&($2 in v)&&$3<'$t'{if($1$2 in a) share[$2]++}END{print length(share),'$total',length(share)/'$total'}' snp.shared.braingvex-rosmap.sig pqtl.braingvex.nominal.chrall pqtl.rosmap.nominal.chrall
done

