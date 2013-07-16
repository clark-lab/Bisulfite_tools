#!/bin/bash
source /etc/profile.d/modules.sh

module load gi/samtools/0.1.19
module load gi/bowtie/2.1.0
module load gi/bismark/0.7.12
module load gi/picard-tools/1.91

GENOMES=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation

OUTFOLDER="$1"_"$SGE_TASK_ID"
FW="$1"_R1_"$SGE_TASK_ID".fastq.gz
RV="$1"_R2_"$SGE_TASK_ID".fastq.gz

#map
bismark -p 4 --bowtie2 -X 1000 --unmapped --ambiguous --gzip --bam -o "$OUTFOLDER" "$GENOMES"/"$2"/bismark_2_sorted/ -1 "$FW" -2 "$RV"

#remove temporary fastqs
rm $FW $RV

#reheader bam
java -jar "$PICARD_HOME"/ReorderSam.jar I="$OUTFOLDER"/"$FW"_bismark_bt2_pe.bam O="$OUTFOLDER"_unsorted.bam R="$GENOMES"/"$2"/bismark_2_sorted/"$2".fa
rm "$OUTFOLDER"/"$FW"_bismark_bt2_pe.bam;

#sort
samtools sort "$OUTFOLDER"_unsorted.bam "$OUTFOLDER"_names
rm "$OUTFOLDER"_unsorted.bam

#remove trailing /1 & /2s from read names
samtools view -h "$OUTFOLDER"_names.bam | awk -F "\t" 'BEGIN {OFS="\t"}{gsub("/[12]", "", $1); print $0}' | samtools view -Sb - > "$OUTFOLDER".bam
rm "$OUTFOLDER"_names.bam

rm -rf "$OUTFOLDOR"
