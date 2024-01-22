############## for clump
##### 7 gene' beta isn't sig. 23 in 30 replicated
awk -F "\t" '{print $4,$1,$2,$3}' ../../metaxcan.gene.bed |while read a b c d;do intercept_b=`awk -F "\t" 'NR==2{print $2}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/ld-clump/eQTL-rQTL/chr"$b"/"$a"/slope.txt`;intercept_p=`awk -F "\t" 'NR==2{print $5}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/ld-clump/eQTL-rQTL/chr"$b"/"$a"/slope.txt`;awk -F "\t" 'NR==3{print "'$a'""\t"$2"\t"$3"\t"$5"\t""'$intercept_b'""\t""'$intercept_p'"}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/ld-clump/eQTL-rQTL/chr"$b"/"$a"/slope.txt;done >metaxcan.gene.clump.eQTL-rQTL


########### for susie
awk -F "\t" '{print $4,$1,$2,$3}' ../../metaxcan.gene.bed |while read a b c d;do intercept_b=`awk -F "\t" 'NR==2{print $2}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/fine-mapping/eQTL-rQTL/chr"$b"/"$a"/slope.txt`;intercept_p=`awk -F "\t" 'NR==2{print $5}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/fine-mapping/eQTL-rQTL/chr"$b"/"$a"/slope.txt`;awk -F "\t" 'NR==3{print "'$a'""\t"$2"\t"$3"\t"$5"\t""'$intercept_b'""\t""'$intercept_p'"}' /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/fine-mapping/eQTL-rQTL/chr"$b"/"$a"/slope.txt;done >metaxcan.gene.susie.eQTL-rQTL


####### merge
Rscript merge.r
awk -F "\t" 'BEGIN{OFS="\t"}ARGIND==1{a[$7]=$6}ARGIND==2{if($1 in a || FNR==1){$1=$1"\t"a[$1];print}}' /zs32/data-analysis/liucy_group/jiangyi/database/gencode.v19/gene2coor one-sample-MR-eQTL-rQTL.txt.tmp >one-sample-MR-eQTL-rQTL.txt

########### get F-statistics for weak instrument bias test
#### for clump
grep "Weak instruments"  /vg_sklmgdata_hw_01/data/liangqiuman/psychencode/another_analysis/code/4.MR/one-sample-MR/ld-clump/eQTL-rQTL/chr*/*/diagnostics.txt|awk '{split($1,a,"/");sub(/chr/,"",a[1]);print a[1]"\t"a[2]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' >weak-instrument.diagnosis.ld-clump.txt

sed -i '1ichr\tgene\tdf1\tdf2\tstatistic\tpvalue' weak-instrument.diagnosis.ld-clump.txt

