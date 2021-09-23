
## check risk scores are sorted in bp order and no header

#!/bin/bash
for i in {1..22}
do
  #awk -v chr=$i '($1 == chr)' 23andmePGC_SBayesR_forupdog.txt > scores_chr"$i".txt ## extract the polygenic risk scores for each chromosome
  awk -v chr=$i '($1 == chr)' 23andmePGC_pTclump_forupdog.txt > scores_chr"$i".txt ## extract the polygenic risk scores for each chromosome
  #awk -v chr=$i '($1 == chr)' 23andmePGC_Jun10_GCldsc_withN_2cohort_pst_eff_a1_b0.5_phiauto_genomewide.txt > scores_chr"$i".txt ## extract the polygenic risk scores for each chromosome

  if [ $(cat 23andmePGC_SBayesR_forupdog.txt | wc -l) -gt 100000 ]  ## if genome-wide scoring is available for 100k+ variant then split into 1000 variant clumps else split by genome size
  then
    split -l1000 -d -a 3 --numeric=1 scores_chr"$i".txt chunkscores_chr"$i"_  ## split the risk scores in to clumps of 1000 variants
  else
    awk '($3 < 30000000)' scores_chr"$i".txt > chunkscores_chr"$i"_001
    awk '($3 >= 30000000 && $3 < 60000000)' scores_chr"$i".txt > chunkscores_chr"$i"_002
    awk '($3 >= 60000000 && $3 < 90000000)' scores_chr"$i".txt > chunkscores_chr"$i"_003
    awk '($3 >= 90000000 && $3 < 120000000)' scores_chr"$i".txt > chunkscores_chr"$i"_004
    awk '($3 >= 120000000 && $3 < 150000000)' scores_chr"$i".txt > chunkscores_chr"$i"_005
    awk '($3 >= 150000000 && $3 < 180000000)' scores_chr"$i".txt > chunkscores_chr"$i"_006
    awk '($3 >= 180000000 && $3 < 210000000)' scores_chr"$i".txt > chunkscores_chr"$i"_007
    awk '($3 >= 210000000 && $3 < 240000000)' scores_chr"$i".txt > chunkscores_chr"$i"_008
    awk '($3 >= 240000000)' scores_chr"$i".txt > chunkscores_chr"$i"_009
    find . -name "chunkscores_chr"$i"_*" -empty -delete ## removes empty chunks
  fi

  #qsub -t 2-2 -l h_rt=4:00:00 -N updog_chr"$i" ./updog.sh "$i"  ## test run of the first chunk
  qsub -t 1-$(ls chunkscores_chr"$i"_* | wc -l) -l h_rt=4:00:00 -N updog_chr"$i" ./updog.sh "$i" ## submit each clump on each chromosome as a job running ./updog.sh
done
