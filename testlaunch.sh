#!/bin/bash

    echo "   _   _          _ "
    echo "  | | | |        | | "
    echo "  | | | |____  _ | | ___  ___ "
    echo "  | | | |  _ \/ || |/ _ \/ _ \ "
    echo "  | |_| | | )( (_| ( (_|( (_| | "
    echo "  \_____| ||_/\____|\___/\__  | "
    echo "        | |              ___| | "
    echo "        |_|     v1.0    (_____| "
    echo ""

while getopts ":ht:l:s:r:" opt; do
  case ${opt} in
    h)
    echo "Use [-t] for location of bed/bim/fam format test data"
    echo "Use [-l] for location of bed/bim/fam format ld reference for summary statistics"
    echo "Use [-s] for location of summary statistics"
    echo "Use [-r] for location of risk scores"
    exit 0
      ;;
    t) # process option t
      testloc=$OPTARG
      ;;
    l) # process option l
      ldloc=$OPTARG
      ;;
    s) # process option s
      sumstats=$OPTARG
      ;;
    r) # process option r
      scores=$OPTARG
      ;;
    \?) 
      echo "Invalid Option: -$OPTARG" 1>&2
      echo ""
      echo "Use [-t] for location of bed/bim/fam format test data"
      echo "Use [-l] for location of bed/bim/fam format ld reference for summary statistics"
      echo "Use [-s] for location of summary statistics"
      echo "Use [-r] for location of risk scores"
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))


if [ -z "${testloc}" ] || [ -z "${ldloc}" ] || [ -z "${sumstats}" ] || [ -z "${scores}" ]; then
  echo "Missing file location(s)"
  if [ -z "${testloc}" ]; then
      echo "Use [-t] for location of bed/bim/fam format test data"
  fi
  if [ -z "${ldloc}" ]; then
      echo "Use [-l] for location of bed/bim/fam format ld reference for summary statistics"
  fi
  if [ -z "${sumstats}" ]; then
      echo "Use [-s] for location of summary statistics"
  fi
  if [ -z "${scores}" ]; then
      echo "Use [-r] for location of risk scores"
  fi
exit 1
fi

echo $testloc
echo $ldloc
echo $sumstats
echo $scores
