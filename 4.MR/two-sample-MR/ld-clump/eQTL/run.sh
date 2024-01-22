chr=$1
gene=$2
mkdir chr"$chr"
cd chr"$chr"
mkdir $gene
cd $gene
###
cp ../../mr.r .
Rscript mr.r --chr $chr --gene $gene 1>"$x".result 2>"$x".mr.log
