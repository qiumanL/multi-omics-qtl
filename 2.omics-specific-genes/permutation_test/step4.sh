num=$1
### rna cond
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle"$num"/rna
cp /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/scripts/r2.r .
Rscript r2.r rna
### ribo cond
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle"$num"/ribo
cp /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/scripts/r2.r .
Rscript r2.r ribo
### protein cond
cd /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/shuffle"$num"/protein
cp /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/permutation_test/scripts/r2.r .
Rscript r2.r protein

