#!/bin/bash -e
if [ ! -e "$1" ]; then
  echo "$1"" does not exist!";
  echo "Usage: `basename $0` {bam file}"
  exit 1;
fi

samtools view -f2 -q30 "$1" | cut -f 9 | grep -v "-" | sort -k1,1n | uniq -c
