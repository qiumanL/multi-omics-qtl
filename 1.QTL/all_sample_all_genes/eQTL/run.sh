
parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf ../genotype/genotypes.all.chr{}.realign.vcf.gz --bed ../expression/log2cpm.fgene.fsample.qn.realign.lmPEER30.norm.sub.chr{}.bed.gz --std-err --out eqtl.nominal.chr{} --nominal 0.05 \> log/log.nominal.all.chr{}  ::: {1..22}
## Permutation pass
parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf ../genotype/genotypes.all.chr{}.realign.vcf.gz --bed ../expression/log2cpm.fgene.fsample.qn.realign.lmPEER30.norm.sub.chr{}.bed.gz --std-err --out eqtl.permute.chr{} --permute 10000 \> log/log.permute.all.chr{}  ::: {1..22}
cat eqtl.permute.chr{1..22} > eqtl.permute
Rscript /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/scripts/qtltools_runFDR_cis.R eqtl.permute 0.05 eqtl.permute
## Conditional analysis
parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf ../genotype/genotypes.all.chr{}.realign.vcf.gz --bed ../expression/log2cpm.fgene.fsample.qn.realign.lmPEER30.norm.sub.chr{}.bed.gz --std-err --out eqtl.conditional.chr{} --mapping eqtl.permute.thresholds.txt --region {} \> log/log.conditional.all.chr{}  ::: {1..22}
## all pairs
parallel -j 6 /zs32/data-analysis/liucy_group/liangqiuman/software/qtltools/bin/QTLtools cis --vcf ../genotype/genotypes.all.chr{}.realign.vcf.gz --bed ../expression/log2cpm.fgene.fsample.qn.realign.lmPEER30.norm.sub.chr{}.bed.gz --std-err --out eqtl.nominal.allpair.chr{} --nominal 1 \> log/log.nominal.allpair.chr{}  ::: {1..22}

###
Rscript qvalue.r
