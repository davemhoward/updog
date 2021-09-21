#!/bin/bash

#$ -cwd
#$ -l h_vmem=16G
# -t 22-22

##SGE_TASK_ID=1
##SGE_TASK_ID=$1

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.0.2
module load igmm/apps/plink/2.03

plink2 \
--bfile ldref/g1000_eur_chr$1 \
--chr $1 \
--from-bp $(head -n1 chunkscores_chr$1_$(printf "%03d" ${SGE_TASK_ID}) | awk '{print $3-250000}') \
--to-bp $(tail -n1 chunkscores_chr$1_$(printf "%03d" ${SGE_TASK_ID}) | awk '{print $3+250000}') \
--make-bed \
--out ~/scratch/g1000_eur_chr$1_$(printf "%03d" ${SGE_TASK_ID})

plink2 \
--pfile ukb/ukb_chr$1 \
--chr $1 \
--from-bp $(head -n1 chunkscores_chr$1_$(printf "%03d" ${SGE_TASK_ID}) | awk '{print $3-250000}') \
--to-bp $(tail -n1 chunkscores_chr$1_$(printf "%03d" ${SGE_TASK_ID}) | awk '{print $3+250000}') \
--exclude ukb/duplicatedSNPs_chr$1.txt \
--make-bed \
--out ukb/ukb_chr$1_$(printf "%03d" ${SGE_TASK_ID})

head -n 1 23andmePGC_sumstats.ma > 23andmePGC_sumstats_chr$1_$(printf "%03d" ${SGE_TASK_ID}).ma
awk 'NR==FNR {FILE1[$2]=$0; next} ($1 in FILE1) {print $0}' \
ukb/ukb_chr$1_$(printf "%03d" ${SGE_TASK_ID}).bim 23andmePGC_sumstats.ma >> 23andmePGC_sumstats_chr$1_$(printf "%03d" ${SGE_TASK_ID}).ma

## chromsome as $1 and task id as chunk
Rscript updog3_for_eddie.R $1 ${SGE_TASK_ID}

rm ~/scratch/g1000_eur_chr$1_$(printf "%03d" ${SGE_TASK_ID})*
rm ukb/ukb_chr$1_$(printf "%03d" ${SGE_TASK_ID})*
rm 23andmePGC_sumstats_chr$1_$(printf "%03d" ${SGE_TASK_ID}).ma
