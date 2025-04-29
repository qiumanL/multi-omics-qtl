### r2设为0.5,p=0.0001
chr=$1
gene=$2
mkdir chr"$chr"
mkdir chr"$chr"/"$gene"
cd chr"$chr"/"$gene"
plink_1.9/plink --bfile genotypes.all.chr"$chr".realign --clump eqtl.allP.forgwas.chr"$chr"."$gene" --clump-kb 1000 --clump-r2 0.5 --clump-snp-field SNP --clump-p1 0.0001 --clump-p2 0.05 --clump-field p

awk 'NR>1{print $3}' plink.clumped >proxy_snp

### r2设为0.5,p=0.001
chr=$1
gene=$2
mkdir chr"$chr"
mkdir chr"$chr"/"$gene"
cd chr"$chr"/"$gene"
plink_1.9/plink --bfile genotypes.all.chr"$chr".realign --clump eqtl.allP.forgwas.chr"$chr"."$gene" --clump-kb 1000 --clump-r2 0.5 --clump-snp-field SNP --clump-p1 0.001 --clump-p2 0.05 --clump-field p

awk 'NR>1{print $3}' plink.clumped >proxy_snp

### r2设为0.5,p=0.01
chr=$1
gene=$2
mkdir chr"$chr"
mkdir chr"$chr"/"$gene"
cd chr"$chr"/"$gene"
plink_1.9/plink --bfile genotypes.all.chr"$chr".realign --clump eqtl.allP.forgwas.chr"$chr"."$gene" --clump-kb 1000 --clump-r2 0.5 --clump-snp-field SNP --clump-p1 0.01 --clump-p2 0.05 --clump-field p

awk 'NR>1{print $3}' plink.clumped >proxy_snp

### r2设为0.01,p=0.05
chr=$1
gene=$2
mkdir chr"$chr"
mkdir chr"$chr"/"$gene"
cd chr"$chr"/"$gene"
plink_1.9/plink --bfile genotypes.all.chr"$chr".realign --clump eqtl.allP.forgwas.chr"$chr"."$gene" --clump-kb 1000 --clump-r2 0.01 --clump-snp-field SNP --clump-p1 0.05 --clump-p2 0.05 --clump-field p

awk 'NR>1{print $3}' plink.clumped >proxy_snp

