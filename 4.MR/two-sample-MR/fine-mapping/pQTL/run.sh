chr=$1
gene=$2
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/fine-mapping/pQTL/chr"$chr"
cd $gene
###
cp ../../mr.pqtl.r .
awk -F "\t" '{split($1,a,"|");print a[2]}' ../"$gene".exp |while read x;do Rscript mr.pqtl.r --chr $chr --gene $gene --protein $x 1>"$x".result 2>"$x".mr.log;done
