#!/bin/bash
echo "        _   _          _ "
echo "       | | | |        | | "
echo "    __@@ |_| |____  _ | | ___  ___ "
echo "   O         |  _ \/ || |/ _ \/ _ \ "
echo "    \_____)  | | )( (_| ( (_|( (_| | "
echo "     U  \___/| ||_/\____|\___/\__  | "
echo "             | |              ___| | "
echo "             |_|     v1.0    (_____| "
echo ""

usage() { 
  echo "  updog requires 4 sources of data:"; \
  echo "  1. Test data. This is the genetic dataset you want to make prediction in to"; \
  echo "  2. Linkage Disequilibrium (ld) reference data matching the ancestry of the"; \
  echo "     summary stats (typically the 1000 genomes or the HRC reference data panels"; \
  echo "  3. The genome-wide summary statistics file"; \
  echo "  4. The genome-wide genetic scores for making predictions (typically created"; \
  echo "     from passing the summary statistics file to packages such as PRS-CS or GCTB,"; \
  echo "     or by using P-value thresholding and clumping)"; \
  echo ""; \
  echo "  To run updog type ./updog and supply the following flags and arguments."; \
  echo "  -t [location of test data]"; \
  echo "     Data should be split by chromosome. Chromosome number should be the last"; \
  echo "     part of the filename prefix and omitted along with the file type. So if the"; \
  echo "     full filename was UK_biobank_chr1.vcf then supply -t UK_biobank_chr"; \
  echo "  -u [filetype of test data]"; \
  echo "     e.g. -u bfile, -u pfile, or -u vcf"; \
  echo "  -l <location of ld reference data for summary statistics>"; \
  echo "     Data should be split by chromosome, chromosome number should be the last"; \
  echo "     part of the filename prefix and omitted along with the file type. So if the"; \
  echo "     full file name was 1000g_eur_chr1.vcf then type -t 1000g_eur_chr"; \
  echo "  -m <filetype of ld reference for summary statistics>"; \
  echo "     e.g. -m bfile, -m pfile, or -m vcf"; \
  echo "  -s <location of summary statistics> {optional flag -n}"; \
  echo "     File should be genome-wide, space separated, with a header row. The first 4"; \
  echo "     columns must contain SNP Name, A1 allele, A2 allele, Effect Size. The -n"; \
  echo "     flag indicates that effect sizes are odds ratios else beta coefficients are"; \
  echo "     assumed. The header row is ignored along with any additional columns."; \
  echo "  -g <location of genetic scores> {optional flag -i}"; \
  echo "     File should be genome-wide, space separated, with a header row. The first 6"; \
  echo "     columns must contain chromosome number, basepair position, SNP Name, A1"; \
  echo "     allele, A2 allele, Genetic Score. The -i flag indicates that genetic scores"; \
  echo "     are odds ratios else beta coefficients are assumed. The header row is"; \
  echo "     ignored along with any additional columns"; \
  echo "  -p <location of plink environment module>"; \
  echo "     i.e. what you enter after entering module load on the command line"; \
  echo "  -r <location of R environment module>"; \
  echo "     i.e. what you enter after entering module load on the command line"; \
  echo "     Use R version ≥ 3.0.2"; \
  echo "  -o <name for output>"; }

sumstatsOR=0
scoresOR=0

while getopts ":hnit:l:s:r:u:m:o:p:g:" opt; do
  case ${opt} in
    h)
      usage
      exit 0
      ;;
    t) # process option t
      testloc=$OPTARG
      ;;
    u) # process option u
      testtype=$OPTARG
      ;;
    l) # process option l
      ldloc=$OPTARG
      ;;
    m) # process option m
      ldtype=$OPTARG
      ;;
    s) # process option s
      sumstats=$OPTARG
      ;;
    n) # process option n
      sumstatsOR=1
      ;;
    g) # process option g
      scores=$OPTARG
      ;;
    i) # process option i
      scoresOR=1
      ;;
    p) # process option p
      plinkloc=$OPTARG
      ;;
    r) # process option r
      rloc=$OPTARG
      ;;
    o) # process option o
      outname=$OPTARG
      ;;
    \?) 
      echo "Invalid Option: -$OPTARG" 1>&2
      echo ""
      usage
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))


if [ -z "${testloc}" ] || [ -z "${ldloc}" ] || [ -z "${sumstats}" ] || \
 [ -z "${scores}" ] || [ -z "${testtype}" ] || [ -z "${ldtype}" ] || \
 [ -z "${plinkloc}" ] || [ -z "${rloc}" ] || [ -z "${outname}" ]; then
  usage
  echo ""
  echo "  Missing flag(s) or  argument(s):"
  if [ -z "${testloc}" ]; then
    echo "  -t"
  fi
  if [ -z "${testtype}" ]; then
    echo "  -u"
  fi
  if [ -z "${ldloc}" ]; then
    echo "  -l"
  fi
  if [ -z "${ldtype}" ]; then
    echo "  -m"
  fi
  if [ -z "${sumstats}" ]; then
    echo "  -s"
  fi
  if [ -z "${scores}" ]; then
    echo "  -g"
  fi
  if [ -z "${plinkloc}" ]; then
    echo "  -p"
  fi
  if [ -z "${rloc}" ]; then
    echo "  -r"
  fi
  if [ -z "${outname}" ]; then
    echo "  -o"; \
  fi
  exit 1
fi

if [ -e $(eval echo genomewidescores_${outname}.txt) ]; then
  echo ""
  echo "  Output already found for genomewidescores_${outname}.txt"
  echo "  Either delete, move or rename genomewidescores_${outname}.txt and resubmit."
  echo ""
  exit 1
fi

echo "  Examining test data:"; \
echo "  ${testloc}*${testtype}"; \
if [ ${testtype} == "bfile" ]; then
  if [ ! -e $(eval echo ${testloc}22.bed) ]; then
    echo "  Test data not found when looking for ${testloc}22.bed"
    exit 1
  fi
  if [ ! -e $(eval echo ${testloc}22.bim) ]; then
    echo "  Test data not found when looking for ${testloc}22.bim"
    exit 1
  fi
  if [ ! -e $(eval echo ${testloc}22.fam) ]; then
    echo "  Test data not found when looking for ${testloc}22.fam"
    exit 1
  fi
  echo "  Test data found with $(cat ${testloc}22.fam | wc -l) individuals and $(cat ${testloc}*.bim | wc -l) variants"; \
fi
if [ ${testtype} == "pfile" ]; then
  if [ ! -e $(eval echo ${testloc}22.pgen) ]; then
    echo "  Test data not found when looking for ${testloc}22.pgen"
    exit 1
  fi
  if [ ! -e $(eval echo ${testloc}22.pvar) ]; then
    echo "  Test data not found when looking for ${testloc}22.pvar"
    exit 1
  fi
  if [ ! -e $(eval echo ${testloc}22.psam) ]; then
    echo "  Test data not found when looking for ${testloc}22.psam"
    exit 1
  fi
  echo "  Test data found with $(cat ${testloc}22.psam | wc -l) individuals and $(cat ${testloc}*.pvar | wc -l) variants"; \
fi
if [ ${testtype} == "vcf" ]; then
  if [ ! -e $(eval echo ${testloc}22.vcf) ]; then
    echo "  Test data not found when looking for ${testloc}22.vcf"
    exit 1
  fi
fi

echo ""
echo "  Examining ld reference data:"; \
echo "  ${ldloc}*${ldtype}"; \
if [ ${ldtype} == "bfile" ]; then
  if [ ! -e $(eval echo ${ldloc}22.bed) ]; then
    echo "  Reference ld data not found when looking for ${ldloc}22.bed"
    exit 1
  fi
  if [ ! -e $(eval echo ${ldloc}22.bim) ]; then
    echo "  Reference ld data not found when looking for ${ldloc}22.bim"
    exit 1
  fi
  if [ ! -e $(eval echo ${ldloc}22.fam) ]; then
    echo "  Reference ld data not found when looking for ${ldloc}22.fam"
    exit 1
  fi
  echo "  LD reference data found with $(cat ${ldloc}22.fam | wc -l) individuals and $(cat ${ldloc}*.bim | wc -l) variants"; \
fi
if [ ${ldtype} == "pfile" ]; then
  if [ ! -e $(eval echo ${ldloc}22.pgen) ]; then
    echo "  Reference ld data not found when looking for ${ldloc}22.pgen"
    exit 1
  fi
  if [ ! -e $(eval echo ${ldloc}22.pvar) ]; then
    echo "  Reference ld data not found when looking for ${ldloc}22.pvar"
    exit 1
  fi
  if [ ! -e $(eval echo ${ldloc}22.psam) ]; then
    echo "  Reference ld data not found when looking for ${ldloc}22.psam"
    exit 1
  fi
  echo "  LD reference data found with $(cat ${ldloc}22.psam | wc -l) individuals and $(cat ${ldloc}*.pvar | wc -l) variants"; \
fi
if [ ${ldtype} == "vcf" ]; then
  if [ ! -e $(eval echo ${ldloc}22.vcf) ]; then
    echo "  Reference ld data not found when looking for ${ldloc}22.vcf"
    exit 1
  fi
fi
echo ""

if [ -e $(eval echo ${sumstats}) ]; then
  echo "  Summary statistics found for $(cat ${sumstats} | tail -n +2 | wc -l) variants with a mean effect size of $(awk '{ total += $4 } END { print total/NR }' $sumstats)"
  if [ "$sumstatsOR" -eq 1 ]; then
    echo "  Odds ratio flag detected for effect sizes - converting by taking the natural log"
  fi
else
  echo "  Summary statistics not found when looking for ${sumstats}"
  exit 1
fi
echo ""

if [  -e $(eval echo ${scores}) ]; then
  echo "  Genetic scores found for $(cat ${scores} | tail -n +2 | wc -l) variants with a mean genetic score of $(awk '{ total += $6 } END { print total/NR }' $scores)"
  if [ "$scoresOR" -eq 1 ]; then
    echo "  Odds ratio flag detected for genetic scores - converting by taking the natural log"
  fi
else
  echo "  Genetic scores not found when looking for ${scores}"
  exit 1
fi

testplink=$( { module load ${plinkloc}; } 2>&1 )
if [[ "$testplink" == *"ERROR"* ]]; then
  echo "Environment module not found for plink at ${plinkloc}"
  echo "Try $ module avail to find the plink module"
  exit 1
fi

testr=$( { module load ${rloc}; } 2>&1 )
if [[ "$testr" == *"ERROR"* ]]; then
  echo "Environment module not found for R at ${rloc}"
  echo "Try $ module avail to find the R module (≥ 3.0.2)"
  exit 1
fi

echo ""
echo "  Initial checks complete. Moving on to chunking and job submission."
echo ""

if [ ! -d "temp" ] 
then
  mkdir temp
fi
if [ ! -d "logs" ] 
then
  mkdir logs
fi

for i in {1..22}
do
  outfilescores=temp/scores_${outname}_chr"$i".txt
  outfilechunkscores=temp/chunkscores_${outname}_chr"$i"_
  if [ "$scoresOR" -eq 1 ]; then  ## if scores are OR then convert using log(OR)
  awk -v chr=$i '($1 == chr)' ${scores} | sort -k2 -n | awk '{print $1,$2,$3,$4,$5,log($6)}' > $outfilescores ## extract the polygenic risk scores for each chromosome
  else
  awk -v chr=$i '($1 == chr)' ${scores} | sort -k2 -n > $outfilescores ## extract the polygenic risk scores for each chromosome
  fi
  if [ $(cat ${scores} | wc -l) -gt 100000 ]; then  ## if genome-wide scoring is available for 100k+ variant then split into 1000 variant clumps else split by genome size
    split -l1000 -d -a 3 --numeric=1 $outfilescores $outfilechunkscores ## split the risk scores in to clumps of 1000 variants
  else
    count=0
    if [ $(awk '($2 < 30000000)' $outfilescores | wc -l) -gt 0 ]; then
      count=$((count+1))
      awk '($2 < 30000000)' $outfilescores > "$outfilechunkscores"00"$count"
    fi
    if [ $(awk '($2 >= 30000000 && $2 < 60000000)' $outfilescores | wc -l) -gt 0 ]; then
      count=$((count+1))
      awk '($2 >= 30000000 && $2 < 60000000)' $outfilescores > "$outfilechunkscores"00"$count"
    fi
    if [ $(awk '($2 >= 60000000 && $2 < 90000000)' $outfilescores | wc -l) -gt 0 ]; then
      count=$((count+1))
      awk '($2 >= 60000000 && $2 < 90000000)' $outfilescores > "$outfilechunkscores"00"$count"
    fi
    if [ $(awk '($2 >= 90000000 && $2 < 120000000)' $outfilescores | wc -l) -gt 0 ]; then
      count=$((count+1))
      awk '($2 >= 90000000 && $2 < 120000000)' $outfilescores > "$outfilechunkscores"00"$count"
    fi
    if [ $(awk '($2 >= 120000000 && $2 < 150000000)' $outfilescores | wc -l) -gt 0 ]; then
      count=$((count+1))
      awk '($2 >= 120000000 && $2 < 150000000)' $outfilescores > "$outfilechunkscores"00"$count"
    fi
    if [ $(awk '($2 >= 150000000 && $2 < 180000000)' $outfilescores | wc -l) -gt 0 ]; then
      count=$((count+1))
      awk '($2 >= 150000000 && $2 < 180000000)' $outfilescores > "$outfilechunkscores"00"$count"
    fi
    if [ $(awk '($2 >= 180000000 && $2 < 210000000)' $outfilescores | wc -l) -gt 0 ]; then
      count=$((count+1))
      awk '($2 >= 180000000 && $2 < 210000000)' $outfilescores > "$outfilechunkscores"00"$count"
    fi
    if [ $(awk '($2 >= 210000000 && $2 < 240000000)' $outfilescores | wc -l) -gt 0 ]; then
      count=$((count+1))
      awk '($2 >= 210000000 && $2 < 240000000)' $outfilescores > "$outfilechunkscores"00"$count"
    fi
    if [ $(awk '($2 >= 240000000)' $outfilescores | wc -l) -gt 0 ]; then
      count=$((count+1))
      awk '($2 >= 240000000)' $outfilescores > "$outfilechunkscores"00"$count"
    fi
  fi

  arguments="$i $testloc $testtype $ldloc $ldtype $sumstats $scores $plinkloc $rloc $outname $sumstatsOR"

  if [ $(cat ${outfilescores} | wc -l) -gt 0 ]; then  ## only submit for chromosomes with variants
#   qsub -t 1-1 -l h_rt=0:10:00 -l h_vmem=16G -N updog_chr"$i" -cwd ./updog.sh ${arguments} ## test run of the first chunk
    qsub -t 1-$(ls "$outfilechunkscores"* | wc -l) -l h_rt=4:00:00 -o logs/ -e logs/ -l h_vmem=16G -N updog_chr"$i" -cwd ./updog.sh ${arguments} 
  fi
done

mergearguments="$testloc $testtype $ldloc $ldtype $sumstats $scores $plinkloc $rloc $outname $sumstatsOR"
qsub -l h_rt=4:00:00 -o logs/ -e logs/ -l h_vmem=8G -N merge_chunks -cwd -hold_jid "updog_chr*" ./merge_chunks.sh ${mergearguments}

echo ""
echo "  That's everything submitted. Go and get yourself a coffee and check back in a while."
echo "  There will either be a results file called genomewidescores_${outname}.txt or if any"
echo "  chunks failed to run there will be a file called resubmitjobs_${outname}. If the"
echo "  resubmitjobs file is there, just enter ./resubmitjobs_${outname} on the command line"
echo "  and go and make another coffee."
echo ""
