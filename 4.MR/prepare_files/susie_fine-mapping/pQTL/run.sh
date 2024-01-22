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
cp ../../scripts/pqtl.sh .
bash pqtl.sh $chr $gene 1>pqtl.log 2>pqtl_err.log
cp ../../scripts/susie.pqtl.r .
awk -F "\t" '{split($1,a,"|");print a[2]}' "$gene".exp |while read x;do Rscript susie.pqtl.r "$gene" $x pqtl 1>pqtl.susie."$x".log 2>pqtl_susie_"$x"_err.log;done
cp ../../scripts/tiqu_cs_snp.pqtl.py .
awk -F "\t" '{split($1,a,"|");print a[2]}' "$gene".exp |while read x;do python tiqu_cs_snp.pqtl.py "$gene" $x pqtl 1>pqtl.tiqu."$x".log 2>pqtl_tiqu_"$x"_err.log;done
