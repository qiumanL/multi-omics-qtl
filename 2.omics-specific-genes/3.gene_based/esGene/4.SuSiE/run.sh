chr=$1
gene=$2
tss=$3
tes=$4
mkdir chr"$chr"
cd chr"$chr"
mkdir $gene
cd $gene
###
cp ../../scripts/genotype.sh .
bash genotype.sh $chr $gene $tss $tes 1>genotype.log 2>genotype_err.log

###
cp ../../scripts/eqtl.sh .
bash eqtl.sh $chr $gene 1>eqtl.log 2>eqtl_err.log

###
cp ../../scripts/susie.r .
Rscript susie.r $gene eqtl 1>qtl.susie.log 2>eqtl_susie_err.log
