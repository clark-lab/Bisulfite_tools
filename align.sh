#!/bin/bash
PATH=/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/bin:/home/aarsta/bin:/home/aarsta/src/homer/bin:"$PATH"
BISMARK=/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/src/bismark_v0.7.7
GENOMES=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation
export PYTHONPATH=/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/python-modules/lib64/python

if [ $# -ne 4 ]
then
  echo "Usage: `basename $0` {Project} {FW} {RV} {genome}"
  exit 1
fi

if [ ! -e "$2" ]; then
  echo "$2"" does not exist!";
  echo "Usage: `basename $0` {Project} {FW} {RV} {genome}"
  exit 1;
fi

if [ ! -e "$3" ]; then
  echo "$3"" does not exist!";
  echo "Usage: `basename $0` {Project} {FW} {RV} {genome}"
  exit 1;
fi

if [ ! -e "$GENOMES"/"$4"/bismark_2/"$4".fa ]; then
  echo "$GENOMES"/"$4"/bismark_2/"$4".fa " does not exist!";
  echo "Usage: `basename $0` {Project} {FW} {RV} {genome}"
  exit 1;
fi

if [ -d "$1" ]; then
  echo "$1 directory already exists!"
  exit 1;
fi

echo "Mapping untrimmed reads"
"$BISMARK"/bismark -p 4 --bowtie2 -X 1000 --unmapped --ambiguous "$GENOMES"/"$4"/bismark_2/ -1 "$2" -2 "$3" -o "$1"_untrimmed

sed 's%/[12]\t%\t%' "$1"_untrimmed/"$2"_bismark_bt2_pe.sam | samtools view -Sb - > "$1"_untrimmed/"$1"_untrimmed_raw.bam
samtools sort "$1"_untrimmed/"$1"_untrimmed_raw.bam "$1"_untrimmed/"$1"_untrimmed
rm "$1"_untrimmed/"$1"_untrimmed_raw.bam;

echo "Trimming reads"
/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/src/trim_galore/trim_galore --paired "$2" "$3"

echo "Mapping trimmed reads"
"$BISMARK"/bismark -p 4 --bowtie2 -X 1000 --unmapped --ambiguous "$GENOMES"/"$4"/bismark_2/ -1 "${2%.fastq.gz}"_val_1.fq.gz -2 "${3%.fastq.gz}"_val_2.fq.gz -o "$1"_trimmed
TKCC20130121_E13VA_L001_12_trimmed/TKCC20130121_E13VA_L001_FW_12_val_1.fq.gz_bismark_bt2_pe.sam

sed 's%/[12]\t%\t%' "$1"_trimmed/"${2%.fastq.gz}"_val_1.fq.gz_bismark_bt2_pe.sam | samtools view -Sb - > "$1"_trimmed/"$1"_trimmed_raw.bam
samtools sort "$1"_trimmed/"$1"_trimmed_raw.bam "$1"_trimmed/"$1"_trimmed
rm "$1"_trimmed/"$1"_trimmed_raw.bam;


