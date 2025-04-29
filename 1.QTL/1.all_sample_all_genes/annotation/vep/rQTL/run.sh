PERL5LIB=""
#### prepare vcf
awk -F "\t" 'BEGIN{OFS="\t"}ARGIND==1{a[$2]=$1}ARGIND==2{split($8,x,"_");if($8 in a){print x[1]"\t"x[2]"\t"a[$8]"\t"x[3]"\t"x[4]"\t.\t.\t."}}' snp2rs rqtl.permute.qvalue.sig.txt|sort -t$'\t' -k1n -k2n >rqtl.vcf

####
vep -i rqtl.vcf --fork 4 -o rqtl.out.vcf --assembly GRCh37 --cache --cache_version 105 --dir /database/VEP/hg19/ --offline --fasta /database/VEP/hg19/Homo_sapiens.GRCh37.dna.primary_assembly.fa  --force_overwrite

