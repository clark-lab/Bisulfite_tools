#!/bin/bash -e
TOOLS=`readlink -f "${0%/*}"`

if [ -d "smash" ]; then
  echo "'smash' directory already exists!"
  exit 1;
fi

ln -s /share/Temp/aarsta/Alignment ./smash

basedir="$PWD"

for file in TKCC*;
do
    cd smash;
    "$TOOLS"/align_lane_cluster.sh "$file" "$basedir"/"$file"/"$file"_R1.fastq.gz "$basedir"/"$file"/"$file"_R2.fastq.gz hg19 &> "$file".log &
    cd ..;
done

