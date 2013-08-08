#!/bin/bash -e
source /etc/profile.d/modules.sh

module load kevyin/java/1.7.0_25
module load gi/samtools/0.1.19
module load gi/bowtie/2.1.0
module load gi/bismark/0.8.3
module load gi/picard-tools/1.91
JAVA="java -Djava.io.tmpdir=/share/Temp"

GENOMES=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation

OUTFOLDER="$1"_"$SGE_TASK_ID"
FW="$1"_"$SGE_TASK_ID".fastq.gz

#map
echo `date`" - Beginning alignment of chunk no ""$SGE_TASK_ID" >> ../"$1".alignment.log
bismark -p 4 --bowtie2 --unmapped --ambiguous --gzip --bam -o "$OUTFOLDER" "$GENOMES"/"$2"/bismark_2_sorted/ "$FW"

#remove temporary fastqs
rm $FW

#convert unmapped and ambiguous reads to bam
$JAVA -jar "$PICARD_HOME"/FastqToSam.jar SM="$OUTFOLDER" F1="$OUTFOLDER"/"$FW"_unmapped_reads.txt O="${FW%.fastq.gz}"_unmapped.bam
$JAVA -jar "$PICARD_HOME"/FastqToSam.jar SM="$OUTFOLDER" F1="$OUTFOLDER"/"$FW"_ambiguous_reads.txt O="${FW%.fastq.gz}"_ambiguous.bam

#reheader bam
$JAVA -jar "$PICARD_HOME"/ReorderSam.jar I="$OUTFOLDER"/"$FW"_bismark_bt2.bam O="$OUTFOLDER"_unsorted.bam R="$GENOMES"/"$2"/bismark_2_sorted/"$2".fa
rm "$OUTFOLDER"/"$FW"_bismark_bt2.bam;

#sort
samtools sort "$OUTFOLDER"_unsorted.bam "$OUTFOLDER"_names
rm "$OUTFOLDER"_unsorted.bam

#remove trailing /1 & /2s from read names
samtools view -h "$OUTFOLDER"_names.bam | awk -F "\t" 'BEGIN {OFS="\t"}{gsub("/[12]", "", $1); print $0}' | samtools view -Sb - > "$OUTFOLDER".bam
rm "$OUTFOLDER"_names.bam

rm -rf "$OUTFOLDOR"

echo `date`" - Completed alignment of chunk no ""$SGE_TASK_ID" >> ../"$1".alignment.log
