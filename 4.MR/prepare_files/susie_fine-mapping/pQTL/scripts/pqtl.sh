chr=$1
gene=$2
###
zgrep "$gene" /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/pqtl/pqtl.nominal.allpair.chr"$chr".gz |awk '{print $1"\t"$8"\t"$14"\t"$15}' >"$gene".pqtl
zgrep "$gene" /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/pqtl/expression/sub.expr.chr"$chr".bed.gz |cut -f4,7- >"$gene".exp
awk -F "\t" '{split($1,a,"|");print a[2]}' "$gene".exp |while read x;do awk -F "\t" '{split($1,b,"|");if(b[2]=="'$x'"){print}}' "$gene".pqtl >"$gene"."$x".pqtl;done
awk -F "\t" '{split($1,a,"|");print a[2]}' "$gene".exp |while read x;do awk -F "\t" '{split($1,b,"|");if(b[2]=="'$x'"){print}}' "$gene".exp >"$gene"."$x".exp;done
