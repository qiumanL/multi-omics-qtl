chr=$1
gene=$2
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/fine-mapping/eQTL/chr"$chr"
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/fine-mapping/eQTL/chr"$chr"/$gene
###
cp ../../mr.eqtl.r .
Rscript mr.eqtl.r --chr $chr --gene $gene 1>result 2>mr.log
