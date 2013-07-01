#!/bin/bash -e
Rbin="R --vanilla --quiet --slave"
TOOLS="/share/ClusterShare/software/contrib/Cancer-Epigenetics/Pipelines/Bisulfite_tools"

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
$Rbin -f "$TOOLS"/calculate_depth.R --args "$1" > "$PROJECT".depth
$Rbin -f "$TOOLS"/CpG_bias.R --args "$1" hg19 > "$PROJECT".CpG_bias
"$TOOLS"/fragment_size.sh "$1" > "$PROJECT".fragment

