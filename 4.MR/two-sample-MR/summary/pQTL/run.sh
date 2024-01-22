######## for clump
##### ENSG00000004534 isn't sig
#### gene\tb\tse\pval\tintercept\tpval_intercept
awk -F "\t" '{print $4,$5,$1,$2,$3}' ../../protein.metaxcan.gene.bed |while read a b c d e;do egger_b=`awk -F "\t" 'NR>1{print $4}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/ld-clump/pQTL/"$a"/"$b"/twosampleMR.egger.result.txt`;egger_p=`awk -F "\t" 'NR>1{print $6}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/ld-clump/pQTL/"$a"/"$b"/twosampleMR.egger.result.txt`;awk -F "\t" '$5=="Inverse variance weighted"{print "'$a'""\t""'$b'""\t"$7"\t"$8"\t"$9"\t""'$egger_b'""\t""'$egger_p'"}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/ld-clump/pQTL/"$a"/"$b"/twosampleMR.result.txt;done >pQTL.metaxcan.gene.clump.uvmr


########### for susie
#### 5 gene have pleiotropy, 1 gene' beta isn't sig.
awk -F "\t" '{print $4,$1,$5,$2,$3}' ../../protein.metaxcan.gene.bed |while read a b c d e;do egger_b=`awk -F "\t" 'NR>1{print $4}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/fine-mapping/pQTL/chr"$b"/"$a"/mr/"$c".twosampleMR.egger.result.txt`;egger_p=`awk -F "\t" 'NR>1{print $6}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/fine-mapping/pQTL/chr"$b"/"$a"/mr/"$c".twosampleMR.egger.result.txt`;awk -F "\t" '$5=="Inverse variance weighted"{print "'$a'""\t""'$c'""\t"$7"\t"$8"\t"$9"\t""'$egger_b'""\t""'$egger_p'"}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/fine-mapping/pQTL/chr"$b"/"$a"/mr/"$c".twosampleMR.result.txt;done >pQTL.metaxcan.gene.susie.uvmr


### merge
Rscript merge.r
awk -F "\t" 'BEGIN{OFS="\t"}ARGIND==1{a[$7]=$6}ARGIND==2{if(FNR==1){$1=$1"\tgene name";print}if($1 in a){$1=$1"\t"a[$1];print}}' /zs32/data-analysis/liucy_group/jiangyi/database/gencode.v19/gene2coor two-sample-MR-result.pQTL.txt.tmp >two-sample-MR-result.pQTL.txt

#### sig
Rscript fdr.r
