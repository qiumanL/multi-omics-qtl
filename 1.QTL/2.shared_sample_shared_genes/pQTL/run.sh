## Permutation pass
parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf ../genotype/genotypes.all.chr{}.realign.vcf.gz --bed expression/sub.expr.chr{}.bed.gz --std-err --out pqtl.permute.chr{} --permute 10000 \> log/log.permute.all.chr{}  ::: {1..22}
cat pqtl.permute.chr{1..22} > pqtl.permute
Rscript /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/scripts/qtltools_runFDR_cis.R pqtl.permute 0.05 pqtl.permute
## Conditional analysis
parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf ../genotype/genotypes.all.chr{}.realign.vcf.gz --bed expression/sub.expr.chr{}.bed.gz --std-err --out pqtl.conditional.chr{} --mapping pqtl.permute.thresholds.txt --region {} \> log/log.conditional.all.chr{}  ::: {1..22}
## all pairs
parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf ../genotype/genotypes.all.chr{}.realign.vcf.gz --bed expression/sub.expr.chr{}.bed.gz --std-err --out pqtl.nominal.allpair.chr{} --nominal 1 \> log/log.nominal.allpair.chr{}  ::: {1..22}

###
Rscript qvalue.r
