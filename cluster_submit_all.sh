#!/bin/bash -e
TOOLS=`readlink -f "${0%/*}"`

if [ -d "Alignment" ]; then
  echo "'Alignment' directory already exists!"
  exit 1;
fi

mkdir Alignment

for file in TKCC*;
do
    cd Alignment;
    "$TOOLS"/align_lane_cluster.sh "$file" ../"$file"/"$file"_R1.fastq.gz ../"$file"/"$file"_R2.fastq.gz hg19 &> "$file".log &
    cd ..;
done
