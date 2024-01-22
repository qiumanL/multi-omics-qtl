chr=$1
gene=$2
###
zgrep "$gene" /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/rqtl/rqtl.nominal.allpair.chr"$chr".gz |awk '{print $1"\t"$8"\t"$14"\t"$15}' >"$gene".rqtl
zgrep "$gene" /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/rqtl/expression/sub.expr.chr"$chr".bed.gz |cut -f4,7- >"$gene".exp
