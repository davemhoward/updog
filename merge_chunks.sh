#!/bin/bash

# supplied arguments: $1=$testloc, $2=$testtype, $3=$ldloc, $4=$ldtype
# $5=$sumstats, $6=$scores, $7=$plinkloc, $8=$rloc, $9=$outname, $10=sumstatsOR

. /etc/profile.d/modules.sh
module load "$8"

## pass to merge_chunks.R
Rscript merge_chunks.R ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10}

## if resubmitjobs file created then update to allow execute
if [ -e $(eval echo  resubmitjobs_${9}) ]; then
  chmod 700 resubmitjobs_${9}
fi

## if genomewide scores created then remove chunked files
if [ -e $(eval echo genomewidescores_${9}.txt) ]; then
  rm temp/chunkscores_${9}_chr*
  rm temp/ldref_${9}_chr*
  rm temp/scoreoutput_${9}_chr*
  rm temp/sumstats_${9}_chr*
  rm temp/testdata_${9}_chr*
fi
