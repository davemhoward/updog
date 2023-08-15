#!/bin/bash

# $i $testloc $testtype $ldloc $ldtype $sumstats $scores $plinkloc $rloc $outname $sumstatsOR
# supplied arguments: $i = chr, $2 $testloc, $3=$testtype, $4=$ldloc, $5=$ldtype
# $6=$sumstats, $7=$scores, $8=$plinkloc, $9=$rloc, $10=$outname, $11=sumstatsOR, $12=weighting

module load "$plinkloc"
module load "$rloc"

start=$(head -n1 temp/chunkscores_"$outname"_chr"$i"_$(printf "%03d" ${SLURM_ARRAY_TASK_ID}) | awk '{print $2-500000}')
stop=$(tail -n1 temp/chunkscores_"$outname"_chr"$i"_$(printf "%03d" ${SLURM_ARRAY_TASK_ID}) | awk '{print $2+500000}')

if [ $start -lt 0 ]
then
  start=0
fi

ldarg="--${ldtype} ${ldloc}${i} --chr ${i} --from-bp $start --to-bp $stop --make-bed --out temp/ldref_${outname}_chr${i}_$(printf "%03d" ${SLURM_ARRAY_TASK_ID})"
testarg="--${testtype} ${testloc}${i} --chr ${i} --from-bp $start --to-bp $stop --make-bed --out temp/testdata_${outname}_chr${i}_$(printf "%03d" ${SLURM_ARRAY_TASK_ID})"

echo ""
echo "  Chunking ld reference data"
echo ""
if grep -q "[2]" <<< "$plinkloc"; then ## chunk ldref and test date using plink 2
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

testbim=temp/testdata_${outname}_chr${i}_$(printf "%03d" ${SLURM_ARRAY_TASK_ID}).bim
chunksum=temp/sumstats_${outname}_chr${i}_$(printf "%03d" ${SLURM_ARRAY_TASK_ID}).txt
if [ -f "$testbim" ]; then
  if [ "${sumstatsOR}" -eq 1 ]; then
    awk 'NR==FNR {FILE1[$2]=$0; next} ($1 in FILE1) {print $0}' \
    temp/testdata_${outname}_chr${i}_$(printf "%03d" ${SLURM_ARRAY_TASK_ID}).bim ${sumstats} | awk '{print $1,$2,$3,log($4)}' > ${chunksum}
  else
    awk 'NR==FNR {FILE1[$2]=$0; next} ($1 in FILE1) {print $0}' \
    temp/testdata_${outname}_chr${i}_$(printf "%03d" ${SLURM_ARRAY_TASK_ID}).bim ${sumstats} > ${chunksum}
  fi
fi

## if chunked summary stats are empty then remove testbim
lensumstats=$(wc -l < ${chunksum})
if [ "${lensumstats}" -eq 0 ]; then
  rm "$testbim"
fi

## pass chunk to updog.R if testdata available else remove chunkscore file to enable successful merge
if [ -f "$testbim" ]; then
  Rscript updog.R ${i} ${SLURM_ARRAY_TASK_ID} ${outname} ${weighting}
else
  rm temp/chunkscores_${outname}_chr${i}_$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
fi
