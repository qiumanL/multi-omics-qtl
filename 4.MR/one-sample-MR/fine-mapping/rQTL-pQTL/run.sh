chr=$1
gene=$2
protein=$3
mkdir /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/fine-mapping/rQTL-pQTL/chr"$chr"
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/fine-mapping/rQTL-pQTL/chr"$chr"
mkdir $gene
cd $gene
###
cp ../../one-sample-mr-rqtl-pqtl.r .
Rscript one-sample-mr-rqtl-pqtl.r $chr $gene $protein 1>result 2>mr.log
