chr=$1
gene=$2
protein=$3
mkdir chr"$chr"
cd chr"$chr"
mkdir $gene
cd $gene
mkdir $protein
cd $protein
###
cp ../../../mr.r .
Rscript mr.r --chr $chr --gene $gene --protein $protein 1>"$x".result 2>"$x".mr.log
