chr=$1
gene=$2
###
zgrep "$gene" pqtl.nominal.allpair.chr"$chr".gz |awk '{split($1,a,"|");print a[1]"\t"$8"\t"$14"\t"$15}' >"$gene".pqtl
zgrep "$gene" sub.expr  >"$gene".exp
