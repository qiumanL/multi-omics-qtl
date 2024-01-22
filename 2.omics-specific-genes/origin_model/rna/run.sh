## prepare data
cd data
bash run.sh
cd ..
mkdir log

## build prediction models
## revising the main_dir of scr/model_training.train.r and scr/weights_cov.train.r to current path
parallel -j 10 Rscript scr/model_training.train.r {} \&\> log/model_training.train.{} ::: {1..75}
Rscript scr/weights_cov.train.r 

## predict gene expression to get GReX
python /zs32/home/yjiang/softwares/PrediXcan/Software/PrediXcan.py --predict --dosages data/ --dosages_prefix genotypes.sub.addCM --samples sample.order.fmt --weights output/braingvex.db --output_prefix output/braingvex.predict
cat output/braingvex.predict_predicted_expression.txt |cut -f2- > output/braingvex.predict_predicted_expression.txt2

