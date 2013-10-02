#!/bin/bash -e
Rbin="R --vanilla --quiet --slave"
TOOLS=`readlink -f "${0%/*}"`

if [ $# -ne 2 ]
then
  echo "Usage: `basename $0` {bamfile} {genome}"
  exit 1
fi

if [ ! -e "$1" ]; then
  echo "$1"" does not exist!";
  echo "Usage: `basename $0` {bamfile} {genome}"
  exit 1;
fi

PROJECT=${1%.bam}

#metrics
$Rbin -f "$TOOLS"/calculate_depth.R --args "$1" > "$PROJECT".depth
$Rbin -f "$TOOLS"/CpG_bias.R --args "$1" "$2" > "$PROJECT".CpG_bias
"$TOOLS"/fragment_size.sh "$1" > "$PROJECT".fragment

