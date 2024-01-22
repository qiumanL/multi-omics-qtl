### r2设为0.5
chr=$1
gene=$2
mkdir /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/plink_ld-clump/eQTL/chr"$chr"/"$gene"
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/prepare_files/plink_ld-clump/eQTL/chr"$chr"/"$gene"
/zs32/home/yjiang/softwares/plink_1.9/plink --bfile /zs32/data-analysis/liucy_group/jiangyi/database/1000G/1000G_EUR_Phase3_plink/1000G.EUR.QC."$chr" --clump /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/1.qtl/shared/eqtl/for_mr/chr"$chr"/rsid/eqtl.allP.forgwas.chr"$chr"."$gene".rsid --clump-kb 1000 --clump-r2 0.5 --clump-snp-field SNP --clump-p1 0.05 --clump-p2 0.05
