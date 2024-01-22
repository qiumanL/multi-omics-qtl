chr=$1
gene=$2
cd chr"$chr"
cd $gene
###
cp ../../scripts/mr.rqtl.r .
Rscript mr.rqtl.r --chr $chr --gene $gene 1>result 2>mr.log
