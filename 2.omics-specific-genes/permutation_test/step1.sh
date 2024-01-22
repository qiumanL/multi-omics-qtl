num=$1
mkdir shuffle"$num"
mkdir shuffle"$num"/rna
mkdir shuffle"$num"/rna/data
mkdir shuffle"$num"/rna/data/predicted_exp
mkdir shuffle"$num"/rna/data/predicted_exp/output
cd shuffle"$num"/rna/data/predicted_exp/output
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/rna/output/braingvex.predict_predicted_expression.txt2 .
cp ../../../../../scripts/original/shuffle.py .
python3 shuffle.py $num
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/2.predixcan/shuffle
### ribo orig
mkdir shuffle"$num"/ribo
mkdir shuffle"$num"/ribo/data
mkdir shuffle"$num"/ribo/data/predicted_exp
mkdir shuffle"$num"/ribo/data/predicted_exp/output
cd shuffle"$num"/ribo/data/predicted_exp/output
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/ribo/output/braingvex.predict_predicted_expression.txt2 .
cp ../../../../../scripts/original/shuffle.py .
python3 shuffle.py $num
### protein orig
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/2.predixcan/shuffle
mkdir shuffle"$num"/protein
mkdir shuffle"$num"/protein/data
mkdir shuffle"$num"/protein/data/predicted_exp
mkdir shuffle"$num"/protein/data/predicted_exp/output
cd shuffle"$num"/protein/data/predicted_exp/output
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/protein/output/braingvex.predict_predicted_expression.txt2 .
cp ../../../../../scripts/original/shuffle.py .
python3 shuffle.py $num
### rna cond
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/2.predixcan/shuffle
cd shuffle"$num"/rna/data
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/rna/data/geno* .
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/rna/data/gencode.v18.coding.gene.txt
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/rna/data/sub.expr
ln -s predicted_exp/output/sample.order
cp ../../../scripts/conditional/conditional.lm.R .
Rscript conditional.lm.R rna $num
sed -i 's/^\t/geneid\t/g' sub.expr.resid
# scr
cd ../
cp -r ../../scripts/conditional/rna/scr .
sed -i 's#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/predixcan/conditional/non_divided/shuffle/shuffle5/rna#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/2.predixcan/shuffle/shuffle'$num'/rna#g' scr/model_training.r
# log
mkdir log
parallel -j 20 Rscript scr/model_training.r {} \&\> log/model_training.{} ::: {1..51}
### ribo cond
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/2.predixcan/shuffle
cd shuffle"$num"/ribo/data
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/ribo/data/geno* .
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/ribo/data/gencode.v18.coding.gene.txt
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/ribo/data/sub.expr
ln -s predicted_exp/output/sample.order
cp ../../../scripts/conditional/conditional.lm.R .
Rscript conditional.lm.R ribo $num
sed -i 's/^\t/geneid\t/g' sub.expr.resid
# scr
cd ../
cp -r ../../scripts/conditional/ribo/scr .
sed -i 's#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/predixcan/conditional/non_divided/shuffle/shuffle5/rna#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/2.predixcan/shuffle/shuffle'$num'/ribo#g' scr/model_training.r
# log
mkdir log
parallel -j 20 Rscript scr/model_training.r {} \&\> log/model_training.{} ::: {1..51}
### protein cond
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/2.predixcan/shuffle
cd shuffle"$num"/protein/data
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/protein/data/geno* .
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/protein/data/protein.pos* .
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/protein/data/sub.expr
ln -s predicted_exp/output/sample.order
cp ../../../scripts/conditional/conditional.lm.R .
Rscript conditional.lm.R protein $num
sed -i 's/^\t/geneid\t/g' sub.expr.resid
# scr
cd ../
cp -r ../../scripts/conditional/protein/scr .
sed -i 's#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/predixcan/conditional/non_divided/protein#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/2.predixcan/shuffle/shuffle'$num'/protein#g' scr/model_training.r
# log
mkdir log
parallel -j 20 Rscript scr/model_training.r {} \&\> log/model_training.{} ::: {1..51}
