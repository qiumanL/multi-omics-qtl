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
cp ../../scripts/rqtl.sh .
bash rqtl.sh $chr $gene 1>rqtl.log 2>rqtl_err.log
cp ../../scripts/susie.r .
Rscript susie.r $gene rqtl 1>rqtl.susie.log 2>rqtl_susie_err.log
cp ../../scripts/tiqu_cs_snp.py .
python3 tiqu_cs_snp.py $gene rqtl 1>rqtl.tiqu.log 2>rqtl_tiqu_err.log
