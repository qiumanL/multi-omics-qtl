awk '
    BEGIN{FS="\t";OFS="\t"}
    ARGIND==1{order[$1]=FNR+1}
    ARGIND==2{gene[$1]++}
    ARGIND==3{
        if(FNR==1){
            for(i=2;i<=NF;i++){
                if($i in order){
                    ind[order[$i]]=i
                }
            }
            printf($1)
            for(i=2;i<=length(ind)+1;i++){
                n=ind[i]
                printf("\t"$n)
            }
            printf("\n")
        }else if($1 in gene){
            printf($1)
            for(i=2;i<=length(ind)+1;i++){
                n=ind[i]
                printf("\t"$n)
            }
            printf("\n")
        }
    }
' ../../share.samples ../../share.genes origin.expr > sub.expr

awk 'BEGIN{FS="\t";OFS="\t"}NR==FNR&&$5=="gene"{split($8,s,".");a[s[1]]=$1"\t"$2"\t"$3"\t"s[1]"\t"s[1]"\t"$6}NR>FNR&&FNR==1{$1="#chr\tstart\tend\tgene\tgene\tstrand";print $0}NR>FNR&&FNR>1&&($1 in a){$1=a[$1];print $0}' <(zcat /zs32/data-analysis/liucy_group/jiangyi/database/gencode.v19/gencode.v19.annotation.bed.gz) sub.expr |sed -e's/^chr//' > sub.expr.bed  # convert to bed
cat <(awk 'NR==1' sub.expr.bed) <(awk 'NR>1' sub.expr.bed |sort -k1,1n -k2,2n) |bgzip > sub.expr.bed.gz
tabix -p bed sub.expr.bed.gz

## Split by chr
parallel -j 11 tabix sub.expr.bed.gz {} -h \|bgzip \> sub.expr.chr{}.bed.gz ::: {1..22}
for i in {1..22}; do tabix -p bed sub.expr.chr$i.bed.gz; done

