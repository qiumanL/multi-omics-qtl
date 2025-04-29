
## Permutation pass
parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf genotypes.all.chr{}.realign.vcf.gz --bed expression/sub.expr.resid.chr{}.bed.gz --std-err --out pqtl.permute.chr{} --permute 10000 \> log/log.permute.all.chr{}  ::: {1..22}

