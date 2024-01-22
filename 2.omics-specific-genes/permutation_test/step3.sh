num=$1
### rna cond
cd shuffle"$num"/rna
# scr
sed -i 's#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle1/rna#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle'$num'/rna#g' scr/weights_cov.r
Rscript scr/weights_cov.r
cd output
cp ../../../scripts/conditional/writeWeights.py .
python3 writeWeights.py
### get conditional GReX
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle"$num"/rna/data
awk '{print $1"\t"$1}' sample.order >sample.order.fmt
# scr
cd ../
python3 /zs32/home/yjiang/softwares/PrediXcan/Software/PrediXcan.py --predict --dosages data/ --dosages_prefix genotypes.sub.addCM --samples sample.order.fmt --weights output/braingvex.db --output_prefix output/braingvex.predict
cat output/braingvex.predict_predicted_expression.txt |cut -f2- > output/braingvex.predict_predicted_expression.txt2


### ribo cond
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/2.predixcan/shuffle/shuffle"$num"/ribo
# scr
sed -i 's#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle1/ribo#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle'$num'/ribo#g' scr/weights_cov.r
Rscript scr/weights_cov.r
cd output
cp ../../../scripts/conditional/writeWeights.py .
python3 writeWeights.py
### get conditional GReX
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle"$num"/ribo/data
awk '{print $1"\t"$1}' sample.order >sample.order.fmt
# scr
cd ../
python3 /zs32/home/yjiang/softwares/PrediXcan/Software/PrediXcan.py --predict --dosages data/ --dosages_prefix genotypes.sub.addCM --samples sample.order.fmt --weights output/braingvex.db --output_prefix output/braingvex.predict
cat output/braingvex.predict_predicted_expression.txt |cut -f2- > output/braingvex.predict_predicted_expression.txt2


### protein cond
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle"$num"/protein
# scr
sed -i 's#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle1/protein#/vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle'$num'/protein#g' scr/weights_cov.r
Rscript scr/weights_cov.r
cd output
cp ../../../scripts/conditional/writeWeights.py .
python3 writeWeights.py
### protein cond
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle"$num"/protein/data
awk '{print $1"\t"$1}' sample.order >sample.order.fmt
# scr
cd ../
python3 /zs32/home/yjiang/softwares/PrediXcan/Software/PrediXcan.py --predict --dosages data/ --dosages_prefix genotypes.sub.addCM --samples sample.order.fmt --weights output/braingvex.db --output_prefix output/braingvex.predict
cat output/braingvex.predict_predicted_expression.txt |cut -f2- > output/braingvex.predict_predicted_expression.txt2
