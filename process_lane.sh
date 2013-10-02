#!/bin/bash -e
source /etc/profile.d/modules.sh
export PATH=/share/ClusterShare/software/contrib/gi/R/3.0.0/bin/:$PATH

module load kevyin/java/1.7.0_25
module load gi/samtools/0.1.19
module load gi/picard-tools/1.91
module load gi/bismark/0.8.3
module load gi/R/3.0.0
module load fabbus/perl/5.14.2
TOOLS=`readlink -f "${0%/*}"`
BISMARK_OPTIONS="--comprehensive --merge_non_CpG --genome_folder /share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation/"$2"/bismark_2/ --bedgraph --counts --report --gzip --buffer_size 20G"
Rbin="R --vanilla --quiet --slave"
LOGFILE="$1".alignment.log
JAVA="java -Djava.io.tmpdir=/share/Temp"

#merge logs
cat trimmed_split/*/*report.txt >> $LOGFILE

echo `date`" - Merging mapped, unmapped and ambiguous reads" >> $LOGFILE

#merge reads
if [ -e trimmed_split/"$1"_R1_1_unmapped.bam ] #Paired end
then
    BISMARK_OPTIONS="-p --no_overlap --ignore_r2 4 "$BISMARK_OPTIONS
    samtools merge "$1"_R1.unmapped.bam trimmed_split/"$1"_R1_*_unmapped.bam
    samtools merge "$1"_R2.unmapped.bam trimmed_split/"$1"_R2_*_unmapped.bam
    samtools merge "$1"_R1.ambiguous.bam trimmed_split/"$1"_R1_*_ambiguous.bam
    samtools merge "$1"_R2.ambiguous.bam trimmed_split/"$1"_R2_*_ambiguous.bam
else #Single end
    BISMARK_OPTIONS="-s "$BISMARK_OPTIONS
    samtools merge "$1".unmapped.bam trimmed_split/"$1"_*_unmapped.bam
    samtools merge "$1".ambiguous.bam trimmed_split/"$1"_*_ambiguous.bam
fi
rm trimmed_split/*unmapped.bam trimmed_split/*ambiguous.bam
samtools merge - trimmed_split/"$1"_*.bam | $JAVA -jar "$PICARD_HOME"/AddOrReplaceReadGroups.jar I=/dev/stdin O="$1".bam ID="$1" LB="$1" PL=Illumina PU=XXX SM="$1"
rm -rf trimmed_split

echo `date`" - Removing duplicates" >> $LOGFILE
$JAVA -jar "$PICARD_HOME"/MarkDuplicates.jar I="$1".bam O="$1".rmdup.bam M="$1".rmdup.metrics REMOVE_DUPLICATES=TRUE AS=TRUE CREATE_INDEX=TRUE
samtools index "$1".rmdup.bam

echo `date`" - Gathering post alignment statistics" >> $LOGFILE
"$TOOLS"/Bisulfite_stats.sh "$1".rmdup.bam "$2"

#call bismark methylation extractor
echo `date`" - Calling methylation extractor" >> $LOGFILE
samtools sort -n "$1".rmdup.bam "$1".rmdup.name
$(which perl) $(which bismark_methylation_extractor) $BISMARK_OPTIONS "$1".rmdup.name.bam
$Rbin -f "$TOOLS"/bedGraph_to_Rdata.R --args "$1".rmdup.name.bedGraph

