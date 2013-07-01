#!/bin/bash -e
if [ ! -e "$1" ]; then
  echo "$1"" does not exist!";
  echo "Usage: `basename $0` {bam file}"
  exit 1;
fi

samtools sort -n "$1" "${1%.bam}"_name

samtools view -f 0x3 "${1%.bam}"_name.bam | cut -f 9 | sed -e N -e 's/\n/ /' | awk '{if ($1>0) print $1; else if ($2>0) print $2}' | sort -k1,1n | uniq -c
