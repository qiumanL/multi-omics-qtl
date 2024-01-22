cd data #then follow run.sh to prepare data needed

## build prediction models for set
parallel -j 10 Rscript scr/model_training.r {} \&\> log/model_training.{} ::: {1..51}
awk 'FNR>1&&length($NF)!=1{print FILENAME,$0}' weights/*  # check "T.1" allele and revise
Rscript scr/weights_cov.r

python3 /zs32/home/yjiang/softwares/PrediXcan/Software/PrediXcan.py --predict --dosages data/ --dosages_prefix genotypes.sub.addCM.chr --samples sample.order.fmt --weights output/braingvex.db --output_prefix output/braingvex.predict
cat output/braingvex.predict_predicted_expression.txt |cut -f2- > output/braingvex.predict_predicted_expression.txt2

Rscript r2.r rna
