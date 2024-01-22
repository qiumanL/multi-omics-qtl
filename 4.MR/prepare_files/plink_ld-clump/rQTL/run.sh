chr=$1
mkdir chr"$chr"
awk -F "\t" '{print $4}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/bed/share.genes.chr"$chr".bed|while read x;do bash plink.sh "$chr" $x;done
#### 提取proxy snp
awk -F "\t" '{print $4}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/bed/share.genes.chr"$chr".bed|while read x;do awk 'NR>1{print $3}' chr"$chr"/"$x"/plink.clumped >chr"$chr"/"$x"/proxy_snp;done
##### 将 rsid和snpid对应上
awk -F "\t" '{print $4}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/bed/share.genes.chr"$chr".bed|while read x;do awk -F "\t" 'ARGIND==1{a[$1]=$2}ARGIND==2{if($1 in a){print $1"\t"a[$1]}}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/snp2rs/snp2rs.chr"$chr" chr"$chr"/"$x"/proxy_snp >chr"$chr"/"$x"/proxy_snp.snpid;done
