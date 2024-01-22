chr=$1
mkdir chr"$chr"
awk -F "\t" '{print $1,$2}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/pqtl/for_mr/share.genes.protein.chr"$chr".bed|while read x y;do bash plink.sh "$chr" $x $y;done
#### 提取proxy snp
awk -F "\t" '{print $1,$2}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/pqtl/for_mr/share.genes.protein.chr"$chr".bed|while read x y;do awk 'NR>1{print $3}' chr"$chr"/"$x"/"$y"/plink.clumped >chr"$chr"/"$x"/"$y"/proxy_snp;done
##### 将 rsid和snpid对应上
awk -F "\t" '{print $1,$2}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/pqtl/for_mr/share.genes.protein.chr"$chr".bed|while read x y;do awk -F "\t" 'ARGIND==1{a[$1]=$2}ARGIND==2{if($1 in a){print $1"\t"a[$1]}}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/snp2rs/snp2rs.chr"$chr" chr"$chr"/"$x"/"$y"/proxy_snp >chr"$chr"/"$x"/"$y"/proxy_snp.snpid;done
