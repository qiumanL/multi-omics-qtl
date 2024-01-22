num=$1
### rna cond
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test
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
sed -i 's#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle1/rna#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle'$num'/rna#g' scr/model_training.r
# log
mkdir log
parallel -j 20 Rscript scr/model_training.r {} \&\> log/model_training.{} ::: {1..51}
### rna cond
cd shuffle"$num"/rna
# scr
sed -i 's#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle1/rna#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle'$num'/rna#g' scr/weights_cov.r
Rscript scr/weights_cov.r
cd output
cp ../../../scripts/conditional/writeWeights.py .
python3 writeWeights.py

### ribo cond
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/
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
sed -i 's#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle1/ribo#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle'$num'/ribo#g' scr/model_training.r
# log
mkdir log
parallel -j 20 Rscript scr/model_training.r {} \&\> log/model_training.{} ::: {1..51}
### protein cond
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/
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
sed -i 's#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle1/protein#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle'$num'/protein#g' scr/model_training.r
# log
mkdir log
parallel -j 20 Rscript scr/model_training.r {} \&\> log/model_training.{} ::: {1..51}
