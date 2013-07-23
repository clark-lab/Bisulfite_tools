#!/bin/bash
GENOMES=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation
TOOLS=`readlink -f "${0%/*}"`
VERSION="1.0.kevin"

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` {Project}"
  exit 1
fi

PROJECT="$1"
FORWARD="input/""$1""_R1.fastq.gz"
REVERSE="input/""$1""_R2.fastq.gz"
GENOME="hg19"

if [ ! -e "$GENOMES"/"$GENOME"/bismark_2/"$GENOME".fa ]; then
  echo "$GENOMES"/"$GENOME"/bismark_2/"$GENOME".fa " does not exist!";
  echo "Usage: `basename $0` {Project}"
  exit 1;
fi

echo 'Bisulfite Alignment Pipeline v'$VERSION" - spoon edition"
echo 'Starting processing of '$1

echo 'Creating output directories'
mkdir -p output/untrimmed output/trimmed output/trimmed_split

#trim reads
echo 'Preparing raw reads for alignment'
if [[ -z $WAIT_JOB_ID ]]
then #Regular job - the data should be there
  if [ ! -e "$FORWARD" ]; then
    echo "$FORWARD"" does not exist!";
    echo "Usage: `basename $0` {Project}"
    exit 1;
  fi

  if [ ! -e "$REVERSE" ]; then
    echo "$REVERSE"" does not exist!";
    echo "Usage: `basename $0` {Project}"
    exit 1;
  fi

  qsub -v MODULEPATH="$MODULEPATH" -N "$PROJECT""_prep_reads" -wd "$PWD" -e "$PWD"/output -o "$PWD"/output -pe smp 4 -b y "$TOOLS"/prep_reads_spoon.sh "$1" "$FORWARD" "$REVERSE"
else #Spoon job - the data should appear after this hold_jid
  qsub -hold_jid $WAIT_JOB_ID -v MODULEPATH="$MODULEPATH" -N "$PROJECT""_prep_reads" -wd "$PWD" -e "$PWD"/output -o "$PWD"/output -pe smp 4 -b y "$TOOLS"/prep_reads_spoon.sh "$1" "$FORWARD" "$REVERSE"
fi

#alignment script for trimmed
echo 'Submitting mapping jobs to the cluster'
cd output/trimmed_split
for i in `seq 1 20`
do
  FW="$1"_R1_"$i".fastq.gz
  RV="$1"_R2_"$i".fastq.gz
  outdir="$1"_"$i"
  qsub -hold_jid "$1""_prep_reads" -v MODULEPATH="$MODULEPATH" -N "$1"_"$i"_align -pe smp 8 -wd $PWD -b y "$TOOLS"/align.sh $outdir $FW $RV $GENOME
done
cd ..

#cleanup the temporary files
qsub -hold_jid "$1""_*" -N "$1"_cleanup -wd "$PWD"/trimmed_split -b y "$TOOLS"/cleanup_lane_cluster.sh "$1"

#gather stats and call methylation on the aligned reads
qsub -v MODULEPATH="$MODULEPATH" -N "$1"_process_lane -hold_jid "$1""_*" -wd "$PWD" -b y "$TOOLS"/process_lane.sh "$1"

#write summary
if [[ -z $EXEC_JOB_ID ]]
then
  qsub -hold_jid "$1"_process_lane -N "$1"_summarize_lane -m e -M `whoami`@garvan.unsw.edu.au -wd "$PWD" -b y "$TOOLS"/summarise_lane.sh "$1"
else
  qsub -hold_jid "$1"_process_lane -N $EXEC_JOB_ID -m e -M `whoami`@garvan.unsw.edu.au -wd "$PWD" -b y "$TOOLS"/summarise_lane.sh "$1"
fi
