#!/bin/bash -e
Rbin="R-2.16 --vanilla --quiet --slave"
tools="/share/ClusterShare/software/contrib/Cancer-Epigenetics/Pipelines/Bisulfite_tools"
picard="/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/src/picard-tools-1.71"

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` {bamfile}"
  exit 1
fi

if [ ! -e "$1" ]; then
  echo "$1"" does not exist!";
  echo "Usage: `basename $0` {bamfile}"
  exit 1;
fi

PROJECT=${1%.bam}

#metrics
$Rbin -f "$tools"/calculate_depth.R --args "$1" > "$PROJECT".depth
$Rbin -f "$tools"/CpG_bias.R --args "$1" hg19 > "$PROJECT".CpG_bias
"$tools"/fragment_size.sh "$1" > "$PROJECT".fragment

