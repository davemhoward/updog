#!/bin/bash
    echo "                   _ "
    echo "    |\ |\         | | "
    echo "   .. || |____  _ | | ___  ___ "
    echo " o--  \| |  _ \/ || |/ _ \/ _ \ "
    echo "  v__    | | )( (_| ( (_|( (_| | "
    echo "     \___| ||_/\____|\___/\__  | "
    echo "         | |              ___| | "
    echo "         |_|     v1.0    (_____| "
    echo ""
usage() { echo "Use [-t] <location of test data> note: data should be split by chromosome, chromosome number should be the last part of the filename and left off the location along with the file type"; \
      echo "Use [-u] <filetype of test data> e.g. bfile, pfile, vcf, bgen, etc."; \
      echo "Use [-l] <location of ld reference for summary statistics>  note: data should be split by chromosome, chromosome number should be the last part of the filename and left off the location aong with the filetype"; \
      echo "Use [-m] <filetype of ld reference for summary statistics> e.g. bfile, pfile, vcf, bgen, etc."; \
      echo "Use [-s] <location of summary statistics> note: file should be genome-wide and have a header with columns: SNP A1 A2 freq b se p N"; \
      echo "Use [-g] <location of genetic risk scores> note: file should be genoome-wide and have no header with columns: chromosome rsid position A1 A2 beta"; \
      echo "Use [-p] <location of plink environment module> note: this should be what you enter after entering $ module load"; \
      echo "Use [-r] <location of R environment module> note: this should be what you enter after entering $ module load"; \
      echo "Use [-o] <name of output>"; \
      1>&2; }

while getopts ":ht:l:s:r:u:m:o:p:g:" opt; do
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
    g) # process option g
      scores=$OPTARG
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
  echo "Missing argument(s)"
  if [ -z "${testloc}" ]; then
    echo "Use [-t] <location of test data>"
  fi
  if [ -z "${testtype}" ]; then
    echo "Use [-u] <filetype of test data> e.g. bfile, pfile, vcf, bgen, etc."
  fi
  if [ -z "${ldloc}" ]; then
    echo "Use [-l] <location of ld reference for summary statistics>"
  fi
  if [ -z "${ldtype}" ]; then
    echo "Use [-m] <filetype of ld reference for summary statistics> e.g. bfile, pfile, vcf, bgen, etc."
  fi
  if [ -z "${sumstats}" ]; then
    echo "Use [-s] <location of summary statistics> File should have a header with columns: SNP A1 A2 freq b se p N"
  fi
  if [ -z "${scores}" ]; then
    echo "Use [-g] <location of genetic risk scores> File should have no header with columns: chromosome rsid position A1 A2 beta"
  fi
  if [ -z "${plinkloc}" ]; then
    echo "Use [-p] <location of plink environment module> Note: this should be what you enter after entering $ module load"
  fi
  if [ -z "${rloc}" ]; then
    echo "Use [-r] <location of R environment module> Note: this should be what you enter after entering $ module load"
  fi
  if [ -z "${outname}" ]; then
    echo "Use [-o] <name of output>"; \
  fi
  exit 1
fi

echo $testloc
echo $testtype
echo $ldloc
echo $ldtype
echo $sumstats
echo $scores
echo $plinkloc
echo $rloc
echo $outname
