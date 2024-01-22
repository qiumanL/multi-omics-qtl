### r2设为0.5
chr=$1
gene=$2
protein=$3
mkdir /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/plink_ld-clump/pQTL/chr"$chr"/"$gene"
mkdir /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/plink_ld-clump/pQTL/chr"$chr"/"$gene"/"$protein"
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/plink_ld-clump/pQTL/chr"$chr"/"$gene"/"$protein"
/zs32/home/yjiang/softwares/plink_1.9/plink --bfile /zs32/data-analysis/liucy_group/jiangyi/database/1000G/1000G_EUR_Phase3_plink/1000G.EUR.QC."$chr" --clump /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/pqtl/for_mr/chr"$chr"/protein/rsid/pqtl.allP.forgwas.chr"$chr"."$gene"."$protein".rsid --clump-kb 1000 --clump-r2 0.5 --clump-snp-field SNP --clump-p1 0.05 --clump-p2 0.05
