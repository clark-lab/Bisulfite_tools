#!/bin/bash
source /etc/profile.d/modules.sh

module load gi/samtools/0.1.19
module load gi/bowtie/2.1.0
module load gi/bismark/0.7.12
module load gi/picard-tools/1.91

GENOMES=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation

#map
bismark -p 4 --bowtie2 -X 1000 --unmapped --ambiguous --gzip --bam -o "$1" "$GENOMES"/"$4"/bismark_2_sorted/ -1 "$2" -2 "$3"

#reheader bam
java -jar "$PICARD_HOME"/ReorderSam.jar I="$1"/"$2"_bismark_bt2_pe.bam O="$1"_unsorted.bam R="$GENOMES"/"$4"/bismark_2_sorted/"$4".fa
rm "$1"/"$2"_bismark_bt2_pe.bam;

#sort
samtools sort "$1"_unsorted.bam "$1"_names
rm "$1"_unsorted.bam

#remove trailing /1 & /2s from read names
samtools view -h "$1"_names.bam | awk -F "\t" 'BEGIN {OFS="\t"}{gsub("/[12]", "", $1); print $0}' | samtools view -Sb - > "$1".bam
rm "$1"_names.bam

rm -rf "$1"
