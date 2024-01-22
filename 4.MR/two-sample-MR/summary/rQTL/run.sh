######## for clump
awk -F "\t" '{print $4,$1,$2,$3}' ../../ribo.metaxcan.gene.bed |while read a b c d;do egger_b=`awk -F "\t" 'NR>1{print $4}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/ld-clump/rQTL/"$a"/twosampleMR.egger.result.txt`;egger_p=`awk -F "\t" 'NR>1{print $6}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/ld-clump/rQTL/"$a"/twosampleMR.egger.result.txt`;awk -F "\t" '$5=="Inverse variance weighted"{print "'$a'""\t"$7"\t"$8"\t"$9"\t""'$egger_b'""\t""'$egger_p'"}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/ld-clump/rQTL/"$a"/twosampleMR.result.txt;done >rQTL.metaxcan.gene.clump.uvmr

########### for susie
awk -F "\t" '{print $4,$1,$2,$3}' ../../ribo.metaxcan.gene.bed |while read a b c d;do egger_b=`awk -F "\t" 'NR>1{print $4}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/fine-mapping/rQTL/chr"$b"/"$a"/mr/twosampleMR.egger.result.txt`;egger_p=`awk -F "\t" 'NR>1{print $6}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/fine-mapping/rQTL/chr"$b"/"$a"/mr/twosampleMR.egger.result.txt`;awk -F "\t" '$5=="Inverse variance weighted"{print "'$a'""\t"$7"\t"$8"\t"$9"\t""'$egger_b'""\t""'$egger_p'"}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/two-sample-MR/fine-mapping/rQTL/chr"$b"/"$a"/mr/twosampleMR.result.txt;done >rQTL.metaxcan.gene.susie.uvmr

### merge
Rscript merge.r
awk -F "\t" 'BEGIN{OFS="\t"}ARGIND==1{a[$7]=$6}ARGIND==2{if(FNR==1){$1=$1"\tgene name";print}if($1 in a){$1=$1"\t"a[$1];print}}' /zs32/data-analysis/liucy_group/jiangyi/database/gencode.v19/gene2coor two-sample-MR-result.rQTL.txt.tmp >two-sample-MR-result.rQTL.txt
#### sig
Rscript fdr.r


