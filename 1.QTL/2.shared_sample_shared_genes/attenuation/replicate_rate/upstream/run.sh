rm result.pqtl-rqtl result.pqtl-eqtl result.rqtl-eqtl result

IFS=" "
qtls=(pQTL rQTL eQTL)
nqtls=${#qtls[@]}

# Proportion of non-null hypothesis in qtl_2 pair
a=0
while [ $a -lt $nqtls ]; do
    let b=$a+1
    while [ $b -lt $nqtls ]; do
        qtl_1=${qtls[$a]}
        qtl_2=${qtls[$b]}
        echo ${qtl_1}-${qtl_2}
        awk 'ARGIND==1{a[$1]++}ARGIND==2{if($1 in a){print $20}}' <(cat ./${qtl_1}/${qtl_1}.permute.qvalue.sig.txt |sed -e's/|[^\t]*\t/\t/') <(cat ./${qtl_2}/${qtl_2}.permute.qvalue.txt |sed -e's/|[^\t]*\t/\t/') > pvalList.permute.${qtl_1}-${qtl_2} # don't use $qtl.nominal.chr1-22.uniqQTL, because it contains too many significant p values.
        for t in 0.2 0.19 0.18 0.17 0.16 0.15 0.14 0.13 0.12 0.11 0.1 0.08 0.06 0.05 0.03 0.02 0.01 0.008 0.005 0.001; do
            shared=`wc -l pvalList.permute.${qtl_1}-${qtl_2}.qvalue|awk '{print $1}'`;
            x=`awk -F "\t" '$1+0<'$t'{print}' pvalList.permute.${qtl_1}-${qtl_2}.qvalue|wc -l`;
            ratio=$(printf "%.2f" `echo "scale=4;$x/$shared*100"|bc`);
            echo -e "${qtl_1}\t${qtl_2}\t$t\t$x\t$shared\t$ratio" >>result.${qtl_1}-${qtl_2}
        done
        let b=$b+1
    done
    let a=$a+1
done

cat result.pqtl-rqtl result.pqtl-eqtl result.rqtl-eqtl >result
