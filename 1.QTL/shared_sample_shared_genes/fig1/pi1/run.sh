##### upstream
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
        for t in 0.1; do
            awk -F "[\t ]+" 'ARGIND==1&&$22<'$t'{a[$1]=1}ARGIND==2&&($1 in a){print $20}' <(cat ./${qtl_1}/${qtl_1}.permute.qvalue.sig.txt |sed -e's/|[^\t]*\t/\t/') <(cat ./${qtl_2}/${qtl_2}.permute.qvalue.txt |sed -e's/|[^\t]*\t/\t/') > pvalList.${qtl_1}-${qtl_2}.$t # don't use $qtl.nominal.chr1-22.uniqQTL, because it contains too many significant p values.
            Rscript qvalue.R pvalList.${qtl_1}-${qtl_2}.$t 0.01
        done
        let b=$b+1
    done
    let a=$a+1
done

##### downstream
IFS=" "
qtls=(eQTL rQTL pQTL)
nqtls=${#qtls[@]}

## Proportion of non-null hypothesis in qtl_2 pair
a=0
while [ $a -lt $nqtls ]; do
    let b=$a+1
    while [ $b -lt $nqtls ]; do
        qtl_1=${qtls[$a]}
        qtl_2=${qtls[$b]}
        echo ${qtl_1}-${qtl_2}
        for t in 0.1; do
            awk -F "[\t ]+" 'ARGIND==1&&$22<'$t'{a[$1]=1}ARGIND==2&&($1 in a){print $20}' <(cat ./${qtl_1}/${qtl_1}.permute.qvalue.sig.txt |sed -e's/|[^\t]*\t/\t/') <(cat ./${qtl_2}/${qtl_2}.permute.qvalue.txt |sed -e's/|[^\t]*\t/\t/') > pvalList.${qtl_1}-${qtl_2}.$t # don't use $qtl.nominal.chr1-22.uniqQTL, because it contains too many significant p values.
            Rscript qvalue.R pvalList.${qtl_1}-${qtl_2}.$t 0.01
        done
        let b=$b+1
    done
    let a=$a+1
done
