#!/bin/bash

# supplied arguments: $1 = chr, $2=$testloc, $3=$testtype, $4=$ldloc, $5=$ldtype
# $6=$sumstats, $7=$scores, $8=$plinkloc, $9=$rloc, $10=$outname, $11=sumstatsOR

module load "$8"
module load "$9"

start=$(head -n1 temp/chunkscores_${10}_chr"$1"_$(printf "%03d" ${SGE_TASK_ID}) | awk '{print $2-250000}')
stop=$(tail -n1 temp/chunkscores_${10}_chr"$1"_$(printf "%03d" ${SGE_TASK_ID}) | awk '{print $2+250000}')

ldarg="--${5} ${4}${1} --chr ${1} --from-bp $start --to-bp $stop --make-bed --out temp/ldref_${10}_chr${1}_$(printf "%03d" ${SGE_TASK_ID})"
testarg="--${3} ${2}${1} --chr ${1} --from-bp $start --to-bp $stop --make-bed --out temp/testdata_${10}_chr${1}_$(printf "%03d" ${SGE_TASK_ID})"

echo ""
echo "  Chunking ld reference data"
echo ""
if grep -q "plink/2" <<< "$8"; then ## chunk ldref and test date using plink 2
  plink2 $ldarg
  echo ""
  echo "  Chunking test data"
  echo ""
  plink2 $testarg
else ## chunk ldref and test date using plink
  plink $ldarg
  echo ""
  echo "  Chunking test data"
  echo ""
  plink $testarg
fi
echo ""

## chunk summary stats that are in testdata
## if sumstats contain OR then convert using natural log
if [ "${11}" -eq 1 ]; then
awk 'NR==FNR {FILE1[$2]=$0; next} ($1 in FILE1) {print $0}' \
temp/testdata_${10}_chr${1}_$(printf "%03d" ${SGE_TASK_ID}).bim ${6} | awk '{print $1,$2,$3,log($4)}' > temp/sumstats_${10}_chr$1_$(printf "%03d" ${SGE_TASK_ID}).txt
else
awk 'NR==FNR {FILE1[$2]=$0; next} ($1 in FILE1) {print $0}' \
temp/testdata_${10}_chr${1}_$(printf "%03d" ${SGE_TASK_ID}).bim ${6} > temp/sumstats_${10}_chr$1_$(printf "%03d" ${SGE_TASK_ID}).txt
fi

## pass chunk to updog.R
Rscript updog.R ${1} ${SGE_TASK_ID} ${10}
