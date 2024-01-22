chr=$1
gene=$2
mkdir /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/fine-mapping/eQTL-rQTL/chr"$chr"
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/fine-mapping/eQTL-rQTL/chr"$chr"
mkdir $gene
cd $gene
###
cp ../../one-sample-mr-eqtl-rqtl.r .
Rscript one-sample-mr-eqtl-rqtl.r $chr $gene 1>result 2>mr.log
