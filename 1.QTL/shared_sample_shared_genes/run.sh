## Shared samples and genes   183 samples,7822 genes 
head -1 {e,r,p}QTL/expression/origin.expr |awk 'NR%3==2' |cut -f2- |sed -e's/\t/\n/g' |awk '{a[$1]++}END{for(i in a){if(a[i]==3)print i}}' > share.samples
awk 'BEGIN{FS="[\t|]"}FNR>1&&!a[ARGIND,$1]++{print $1}' {e,r,p}QTL/expression/origin.expr |awk '{a[$1]++}END{for(i in a){if(a[i]==3)print i}}' > share.genes

## Check whether the order of sample list matched (matched.)
awk -F "\t" '$1!=$2||$2!=$3||$3!=$4{print}' all.header


### QTL
for i in eQTL rQTL pQTL; do
## nominal qtl
    parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf /zs32/data-analysis/liucy_group/jiangyi/psychENCODE/qtl.subSamples.norm2/genotypes/genotypes.all.chr{}.realign.vcf.gz --bed /zs32/data-analysis/liucy_group/jiangyi/psychENCODE/qtl.subSamples.norm2/"$i"/expression/sub.expr.chr{}.bed.gz --std-err --out "$i"/"$i".nominal.chr{} --nominal 0.05 \> log/"$i"/log.nominal.all.chr{}  ::: {1..22}
### Permutation pass
    parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf /zs32/data-analysis/liucy_group/jiangyi/psychENCODE/qtl.subSamples.norm2/genotypes/genotypes.all.chr{}.realign.vcf.gz --bed /zs32/data-analysis/liucy_group/jiangyi/psychENCODE/qtl.subSamples.norm2/"$i"/expression/sub.expr.chr{}.bed.gz --std-err --out "$i"/"$i".permute.chr{} --permute 10000 \> log/"$i"/log.permute.all.chr{}  ::: {1..22}
    cat "$i"/"$i".permute.chr{1..22} > "$i"/"$i".permute
    Rscript /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/scripts/qtltools_runFDR_cis.R "$i"/"$i".permute 0.05 "$i"/"$i".permute
### Conditional analysis
    parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf /zs32/data-analysis/liucy_group/jiangyi/psychENCODE/qtl.subSamples.norm2/genotypes/genotypes.all.chr{}.realign.vcf.gz --bed /zs32/data-analysis/liucy_group/jiangyi/psychENCODE/qtl.subSamples.norm2/"$i"/expression/sub.expr.chr{}.bed.gz --std-err --out "$i"/"$i".conditional.chr{} --mapping "$i"/"$i".permute.thresholds.txt --region {} \> log/"$i".log.conditional.all.chr{}  ::: {1..22}
### nominal all pairs
    parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf /zs32/data-analysis/liucy_group/jiangyi/psychENCODE/qtl.subSamples.norm2/genotypes/genotypes.all.chr{}.realign.vcf.gz --bed /zs32/data-analysis/liucy_group/jiangyi/psychENCODE/qtl.subSamples.norm2/"$i"/expression/sub.expr.chr{}.bed.gz --std-err --out "$i"/"$i".nominal.allpair.chr{} --nominal 1 \> log/"$i".log.nominal.allpair.chr{}  ::: {1..22}
done
