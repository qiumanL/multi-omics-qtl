Rscript conditional.lm.R protein
sed -i 's/^\t/geneid\t/g' sub.expr.resid

ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/protein/data/genotypes.sub.addCM.* .
ln -s /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/2.omics-specific-genes/origin_model/protein/data/protein.pos* .

## check header
awk -F "\t" 'NR==1{for(i=2;i<=NF;i++){print $i}}' sub.expr.resid >sample.order
paste -d'\t' sample.order /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/2.predixcan/raw/ribo/data/sample.order >sample.order.fmt
awk -F "\t" '$1!=$2{print}' sample.order.fmt 

