
parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf ../genotype/genotypes.all.chr{}.realign.vcf.gz --bed expression/sub.expr.chr{}.bed.gz --std-err --out rqtl.nominal.chr{} --nominal 0.05 \> log/log.nominal.all.chr{}  ::: {1..22}
## Permutation pass
parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf ../genotype/genotypes.all.chr{}.realign.vcf.gz --bed expression/sub.expr.chr{}.bed.gz --std-err --out rqtl.permute.chr{} --permute 10000 \> log/log.permute.all.chr{}  ::: {1..22}
cat rqtl.permute.chr{1..22} > rqtl.permute
Rscript /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/scripts/qtltools_runFDR_cis.R rqtl.permute 0.05 rqtl.permute
## Conditional analysis
parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf ../genotype/genotypes.all.chr{}.realign.vcf.gz --bed expression/sub.expr.chr{}.bed.gz --std-err --out rqtl.conditional.chr{} --mapping rqtl.permute.thresholds.txt --region {} \> log/log.conditional.all.chr{}  ::: {1..22}
## all pairs
parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf ../genotype/genotypes.all.chr{}.realign.vcf.gz --bed expression/sub.expr.chr{}.bed.gz --std-err --out rqtl.nominal.allpair.chr{} --nominal 1 \> log/log.nominal.allpair.chr{}  ::: {1..22}

###
Rscript qvalue.r
