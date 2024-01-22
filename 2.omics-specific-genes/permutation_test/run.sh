i=$1
## 1. shuffle the sampleID of GReX of other two omics
bash step1.sh $i
## 2. follow the steps of conditional_model after shuffle
bash step2.sh $i
## 3. get shuffled GRex based on .db and .cov
bash step3.sh $i
## 4. calculate the R2 between conditional model and 计算 shuffle前后的GReX的相关性平方
bash step4.sh $i
## 5. merge r2 of all 50 shuffles
Rscript scripts/merge.r rna
Rscript scripts/merge.r ribo
Rscript scripts/merge.r protein
## 6. aggregate to the data format for subsequent analysis
Rscript scripts/transfer_table.r rna
awk -F "\t" 'BEGIN{OFS="\t"}ARGIND==1{a[$1]=$2}ARGIND==2{if($1 in a){$1=$1"\t"a[$1];print}}' ../conditional/rna/r2.rna.txt r2.rna.tmp >r2.rna.txt
Rscript scripts/transfer_table.r ribo
awk -F "\t" 'BEGIN{OFS="\t"}ARGIND==1{a[$1]=$2}ARGIND==2{if($1 in a){$1=$1"\t"a[$1];print}}' ../conditional/ribo/r2.ribo.txt r2.ribo.tmp >r2.ribo.txt
Rscript scripts/transfer_table.r protein
awk -F "\t" 'BEGIN{OFS="\t"}ARGIND==1{a[$1]=$2}ARGIND==2{if($1 in a){$1=$1"\t"a[$1];print}}' ../conditional/protein/r2.protein.txt r2.protein.tmp >r2.protein.txt
## 7. revise r2 to observed in header
sed -i 's/r2/observed/g' r2.rna.txt
sed -i 's/r2/observed/g' r2.ribo.txt
sed -i 's/r2/observed/g' r2.protein.txt
