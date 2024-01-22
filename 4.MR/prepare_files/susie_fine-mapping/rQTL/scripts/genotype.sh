chr=$1
gene=$2
tss=$3
tes=$4
###
chrlen=`awk -F "\t" '$1=="chr'$chr'"{print $2}' /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta.fai`
echo $chrlen
if [ $tss -lt 1000000 ];then
	upstream=0
else
	upstream=`expr $tss - 1000000`
fi
echo $upstream
downstream=`expr $tes + 1000000`
echo $downstream
if [ $downstream -gt $chrlen ];then
	downstream=$chrlen
else
	downstream=$downstream
fi
echo $downstream
bcftools filter /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/genotype/genotypes.all.chr"$chr".realign.vcf.gz --regions "$chr":"$upstream"-"$downstream" >"$gene".vcf
vcftools --vcf "$gene".vcf --012 --out "$gene".dosage.vcf

awk '{$1=null;print $0}' "$gene".dosage.vcf.012 > tmp.txt
paste -d" " "$gene".dosage.vcf.012.indv tmp.txt | sed 's/  / /g' > "$gene".dosage.vcf_indv.txt
echo "ID" > snpid.txt
grep -v "^#" "$gene".vcf | cut -f3 >>snpid.txt

cat snpid.txt| tr "\n" " " >snpid2.txt
sed -i 's/[ ]$/\n/g' snpid2.txt
#echo "" >> snpid2.txt
cat snpid2.txt "$gene".dosage.vcf_indv.txt> "$gene".genotype

