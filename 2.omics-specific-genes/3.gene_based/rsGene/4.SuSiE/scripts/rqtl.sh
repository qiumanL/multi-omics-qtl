chr=$1
gene=$2
###
zgrep "$gene" rqtl.nominal.allpair.chr"$chr".gz |awk '{split($1,a,"|");print a[1]"\t"$8"\t"$14"\t"$15}' >"$gene".rqtl
zgrep "$gene" sub.expr  >"$gene".exp
