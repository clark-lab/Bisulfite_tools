#!/bin/bash -e
PATH=/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/bin:/home/aarsta/bin:/home/aarsta/src/homer/bin:"$PATH"
picard="/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/src/picard-tools-1.71"
TOOLS=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Pipelines/Bisulfite_tools/
BISMARK=/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/src/bismark_v0.7.7
BISMARK_OPTIONS="-p --no_overlap --comprehensive --merge_non_CpG --genome_folder /share/ClusterShare/software/contrib/Cancer-Epigenetics//Annotation/hg19/bismark_2/ --bedgraph --counts --report"

##untrimmed data
#samtools merge "$1"_untrimmed.bam */*_untrimmed.bam
#java -jar "$picard"/MarkDuplicates.jar I="$1"_untrimmed.bam O="$1"_untrimmed.rmdup.bam M="$1"_untrimmed.rmdup.metrics REMOVE_DUPLICATES=TRUE AS=TRUE CREATE_INDEX=TRUE
#"$TOOLS"/Bisulfite_stats.sh "$1"_untrimmed.rmdup.bam

#trimmed data
samtools merge "$1"_trimmed.bam */*_trimmed.bam
java -jar "$picard"/MarkDuplicates.jar I="$1"_trimmed.bam O="$1"_trimmed.rmdup.bam M="$1"_trimmed.rmdup.metrics REMOVE_DUPLICATES=TRUE AS=TRUE CREATE_INDEX=TRUE
"$TOOLS"/Bisulfite_stats.sh "$1"_trimmed.rmdup.bam

#call bismark methylation extractor
samtools index "$1"_trimmed.rmdup.bam
samtools view "$1"_trimmed.rmdup.bam > "$1"_trimmed.rmdup.sam
"$BISMARK"/bismark_methylation_extractor $BISMARK_OPTIONS "$1"_trimmed.rmdup.sam

