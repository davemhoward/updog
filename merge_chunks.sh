#!/bin/bash

# supplied arguments: $1=$testloc, $2=$testtype, $3=$ldloc, $4=$ldtype
# $5=$sumstats, $6=$scores, $7=$plinkloc, $8=$rloc, $9=$outname, $10=sumstatsOR $11=jobargs

module load "$rloc"

## pass to merge_chunks.R
Rscript merge_chunks.R ${testloc} ${testtype} ${ldloc} ${ldtype} ${sumstats} ${scores} ${plinkloc} ${rloc} ${outname} ${sumstatsOR} ${jobargs}

## if resubmitjobs file created then update to allow execute
#if [ -e $(eval echo  resubmitjobs_${outname}.qsub) ]; then
#  chmod 700 resubmitjobs_${outname}*
#fi

FILE=resubmitjobs_${outname}.qsub
if test -f "$FILE"; then
  chmod 700 resubmitjobs_${outname}*
fi


## if genomewide scores created then remove chunked files
if [ -e $(eval echo genomewidescores_${outname}.txt) ]; then
  rm temp/chunkscores_${outname}_chr*
  rm temp/ldref_${outname}_chr*
  rm temp/scoreoutput_${outname}_chr*
  rm temp/sumstats_${outname}_chr*
  rm temp/testdata_${outname}_chr*
fi
