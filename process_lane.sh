#!/bin/bash -e
source /etc/profile.d/modules.sh

module load gi/samtools/0.1.19
module load gi/picard-tools/1.91
module load gi/bismark/0.7.12
module load gi/R/3.0.0

TOOLS=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Pipelines/Bisulfite_tools/
BISMARK=/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/src/bismark_v0.7.7
BISMARK_OPTIONS="-p --no_overlap --comprehensive --merge_non_CpG --genome_folder /share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation/hg19/bismark_2_sorted/ --bedgraph --counts --report --gzip --buffer_size 20G"
Rbin="R --vanilla --quiet --slave"

#trimmed data
samtools merge "$1".bam trimmed_split/"$1"_*.bam
java -jar "$PICARD_HOME"/MarkDuplicates.jar I="$1".bam O="$1".rmdup.bam M="$1".rmdup.metrics REMOVE_DUPLICATES=TRUE AS=TRUE CREATE_INDEX=TRUE
"$TOOLS"/Bisulfite_stats.sh "$1".rmdup.bam
rm -rf trimmed_split

#call bismark methylation extractor
samtools index "$1".rmdup.bam
bismark_methylation_extractor $BISMARK_OPTIONS "$1".rmdup.bam
$Rbin -f "$TOOLS"/bedGraph_to_Rdata.R --args "$1".rmdup.bedGraph

#zip up big files
gzip Non_CpG_context_"$f".rmdup.txt
gzip CpG_context_"$f".rmdup.txt

