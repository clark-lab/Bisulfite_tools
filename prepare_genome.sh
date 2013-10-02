#!/bin/bash -e
TOOLS=`readlink -f "${0%/*}"`
Rbin="R --vanilla --quiet --slave"

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` {Genome}"
  exit 1
fi

if [ ! -e "$1".fa ]
then
  echo "$1"".fa does not exist!"
  echo "Usage: `basename $0` {Genome}"
  exit 1;
fi

module load gi/bismark/0.8.3
module load gi/bowtie/2.1.0

BOWTIE2=$(which bowtie2)

bismark_genome_preparation --bowtie2 --path_to_bowtie ${BOWTIE2%bowtie2} ./

module load gi/picard-tools/1.91

java -jar "$PICARD_HOME"/CreateSequenceDictionary.jar R="$1".fa O="$1".dict GENOME_ASSEMBLY="$1"

module load gi/R/3.0.0
$Rbin -f "$TOOLS"/CpG_preparation.R --args "$1"

