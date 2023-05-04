#!/bin/bash

# supplied arguments: $1=$testloc, $2=$testtype, $3=$ldloc, $4=$ldtype, $5=$sumstats, $6=$scores
# $7=$plinkloc, $8=$rloc, $9=$outname, $10=sumstatsOR $11=jobargs, $12=weighting

module load "$rloc"

echo "Chunking complete. Moving on to merging chunks." >> logs/${outname}.txt

## pass to merge_chunks.R
Rscript merge_chunks.R ${testloc} ${testtype} ${ldloc} ${ldtype} ${sumstats} ${scores} ${plinkloc} ${rloc} ${outname} ${sumstatsOR} ${jobargs} ${weighting}

FILE=resubmitjobs_${outname}
if test -f "$FILE"; then
  chmod 700 resubmitjobs_${outname}*
  echo "Some chunks have failed" >> logs/${outname}.txt
  echo "Enter ./resubmitjobs_${outname} on the command line to retry" >> logs/${outname}.txt
fi


## if genomewide scores created then remove chunked files
if [ -e $(eval echo genomewidescores_${outname}.txt) ]; then
  echo "Merging complete. Results written to genomewidescores_${outname}.txt" >> logs/${outname}.txt
  rm temp/chunkscores_${outname}_chr*
  rm temp/ldref_${outname}_chr*
  rm temp/scoreoutput_${outname}_chr*
  rm temp/scores_${outname}_chr*
  rm temp/sumstats_${outname}_chr*
  rm temp/testdata_${outname}_chr*
fi

echo "end time: $(date)" >> logs/${outname}.txt
