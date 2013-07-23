#!/bin/bash
GENOMES=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation
TOOLS=`readlink -f "${0%/*}"`
VERSION="1.0.kevin"

if [ $# -ne 4 ]
then
  echo "Usage: `basename $0` {Project} {FW.fastq.gz} {RV.fastq.gz} {genome}"
  exit 1
fi

if [ ! -e "$2" ]; then
  echo "$2"" does not exist!";
  echo "Usage: `basename $0` {Project} {FW.fastq.gz} {RV.fastq.gz} {genome}"
  exit 1;
fi

if [[ $2 != *.fastq.gz ]]; then
  echo "$2"" filename needs to end in .fastq.gz"; 
  echo "Usage: `basename $0` {Project} {FW.fastq.gz} {RV.fastq.gz} {genome}"
  exit 1;
fi

if [ ! -e "$3" ]; then
  echo "$3"" does not exist!";
  echo "Usage: `basename $0` {Project} {FW.fastq.gz} {RV.fastq.gz} {genome}"
  exit 1;
fi

if [[ $3 != *.fastq.gz ]]; then
  echo "$3"" filename needs to end in .fastq.gz"; 
  echo "Usage: `basename $0` {Project} {FW.fastq.gz} {RV.fastq.gz} {genome}"
  exit 1;
fi

if [ ! -e "$GENOMES"/"$4"/bismark_2/"$4".fa ]; then
  echo "$GENOMES"/"$4"/bismark_2/"$4".fa " does not exist!";
  echo "Usage: `basename $0` {Project} {FW.fastq.gz} {RV.fastq.gz} {genome}"
  exit 1;
fi

if [ -d "$1" ]; then
  echo "$1 directory already exists!"
  exit 1;
fi

echo 'Bisulfite Alignment Pipeline v'$VERSION
echo 'Starting processing of '$1

echo 'Creating output directory'
mkdir "$1" "$1"/untrimmed "$1"/trimmed "$1"/trimmed_split

#trim reads
echo 'Preparing raw reads for alignment'
qsub -v MODULEPATH="$MODULEPATH" -N "$1""_prep_reads" -wd "$PWD" -e "$PWD"/"$1" -o "$PWD"/"$1" -pe smp 4 -b y "$TOOLS"/prep_reads.sh "$1" "$2" "$3"

#alignment script for trimmed
echo 'Submitting mapping jobs to the cluster'
cd "$1"/trimmed_split
for i in `seq 1 20`
do
  FW="$1"_R1_"$i".fastq.gz
  RV="$1"_R2_"$i".fastq.gz
  outdir="$1"_"$i"
  qsub -hold_jid "$1""_prep_reads" -v MODULEPATH="$MODULEPATH" -N "$1"_"$i"_align -pe smp 8 -wd $PWD -b y "$TOOLS"/align.sh $outdir $FW $RV $4
done
cd ..

#cleanup the temporary files
qsub -hold_jid "$1""_*" -N "$1"_cleanup -wd "$PWD"/trimmed_split -b y "$TOOLS"/cleanup_lane_cluster.sh "$1"

#gather stats and call methylation on the aligned reads
qsub -v MODULEPATH="$MODULEPATH" -N "$1"_process_lane -hold_jid "$1""_*" -wd "$PWD" -b y "$TOOLS"/process_lane.sh "$1"

#write summary
qsub -hold_jid "$1"_process_lane -N "$1"_summarize_lane -m e -M `whoami`@garvan.unsw.edu.au -wd "$PWD" -b y "$TOOLS"/summarise_lane.sh "$1"

