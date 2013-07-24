#!/bin/bash -e
source /etc/profile.d/modules.sh

module load kevyin/java/1.7.0_25
module load gi/samtools/0.1.19
module load gi/picard-tools/1.91
module load gi/bismark/0.7.12
module load gi/R/3.0.0

TOOLS=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Pipelines/Bisulfite_tools/
BISMARK=/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/src/bismark_v0.7.7
BISMARK_OPTIONS="-p --no_overlap --comprehensive --merge_non_CpG --genome_folder /share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation/hg19/bismark_2_sorted/ --bedgraph --counts --report --gzip --buffer_size 20G"
Rbin="R --vanilla --quiet --slave"
LOGFILE="$1".alignment.log

#merge logs
cat trimmed_split/*/*_bismark_bt2_PE_report.txt >> $LOGFILE

echo `date`" - Merging mapped, unmapped and ambiguous reads" >> $LOGFILE

#merge reads
samtools merge "$1"_R1.unmapped.bam trimmed_split/"$1"_R1_*_unmapped.bam
samtools merge "$1"_R2.unmapped.bam trimmed_split/"$1"_R2_*_unmapped.bam
samtools merge "$1"_R1.ambiguous.bam trimmed_split/"$1"_R1_*_ambiguous.bam
samtools merge "$1"_R2.ambiguous.bam trimmed_split/"$1"_R2_*_ambiguous.bam
rm trimmed_split/*unmapped.bam trimmed_split/*ambiguous.bam
samtools merge - trimmed_split/"$1"_*.bam | java -jar "$PICARD_HOME"/AddOrReplaceReadGroups.jar I=/dev/stdin O="$1".bam ID="$1" LB="$1" PL=Illumina PU=XXX SM="$1"
rm -rf trimmed_split

echo `date`" - Removing duplicates" >> $LOGFILE
java -jar "$PICARD_HOME"/MarkDuplicates.jar I="$1".bam O="$1".rmdup.bam M="$1".rmdup.metrics REMOVE_DUPLICATES=TRUE AS=TRUE CREATE_INDEX=TRUE
samtools index "$1".rmdup.bam

echo `date`" - Gathering post alignment statistics" >> $LOGFILE
"$TOOLS"/Bisulfite_stats.sh "$1".rmdup.bam

#call bismark methylation extractor
echo `date`" - Calling methylation extractor" >> $LOGFILE
bismark_methylation_extractor $BISMARK_OPTIONS "$1".rmdup.bam
$Rbin -f "$TOOLS"/bedGraph_to_Rdata.R --args "$1".rmdup.bedGraph

#zip up big files
gzip Non_CpG_context_"$f".rmdup.txt
gzip CpG_context_"$f".rmdup.txt

