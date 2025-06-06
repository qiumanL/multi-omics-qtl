
## MESC wiki page: https://github.com/douglasyao/mesc/wiki

## Prepare files
awk -F "\t" 'BEGIN{OFS="\t"}ARGIND==1{split($7,a,".");split($1,chr,"chr");x[a[1]]=chr[2];y[a[1]]=int(($2+$3)/2)}ARGIND==2{if(FNR==1){$1="GENE\tCHR\tGENE_COORD";print $0}else{split($1,pqtl,"|");$1=$1"\t"x[pqtl[1]]"\t"y[pqtl[1]];print}}' gencode.v19.annotation.gene.bed origin.expr |sed -e's/|/__/'  >protein.expr.fmt

## Estimation of expression scores. 
for i in {1..22}; do mesc/run_mesc.py --out expscore/protein.expr.fmt --compute-expscore-indiv --plink-path plink --expression-matrix protein.expr.fmt --exp-bfile bfile/genotypes.all.addCM.chr$i --geno-bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.$i --chr $i --gcta-path gcta64 --tmp tmp; done

## Estimation of h2med. 
mesc/run_mesc.py --h2med 2021PGC3.sczVScontrol.sumstats.gz --exp-chr expscore/protein.expr.fmt --out h2med/protein.expr.fmt.2021PGC3.sczVScontrol
