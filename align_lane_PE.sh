#!/bin/bash -e
QSUB="qsub -q all.q -v MODULEPATH=$MODULEPATH"
GENOMES=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation
TOOLS=`readlink -f "${0%/*}"`
VERSION="1.2.4"

if [ $# -lt 2 ]
then
  echo "Usage: `basename $0` {Project} {Genome} [optional parameters to pass to bismark]"
  exit 1
fi

PROJECT="$1"
LOGFILE="$1".alignment.log
FORWARD="input/""$1""_R1.fastq.gz"
REVERSE="input/""$1""_R2.fastq.gz"
GENOME="$2"

if [ ! -e "$GENOMES"/"$GENOME"/bismark_2/"$GENOME".fa ]; then
  echo "$GENOMES"/"$GENOME"/bismark_2/"$GENOME".fa " does not exist!";
  echo "Usage: `basename $0` {Project}"
  exit 1;
fi

mkdir -p output
echo `date`' - Starting processing of '$1 > output/"$LOGFILE"
echo 'Bisulfite Alignment Pipeline v'$VERSION" - spoon edition" >> output/"$LOGFILE"

echo 'Creating output directories' >> output/"$LOGFILE"
mkdir -p output/untrimmed output/trimmed output/trimmed_split

#trim reads
echo 'Preparing raw reads for alignment' >> output/"$LOGFILE"
if [[ -z $WAIT_JOB_ID ]]
then #Regular job - the data should be there
  if [ ! -e "$FORWARD" ]; then
    echo "$FORWARD"" does not exist!"; >> output/"$LOGFILE"
    echo "Usage: `basename $0` {Project}" >> output/"$LOGFILE"
    exit 1;
  fi

  if [ ! -e "$REVERSE" ]; then
    echo "$REVERSE"" does not exist!"; >> output/"$LOGFILE"
    echo "Usage: `basename $0` {Project}" >> output/"$LOGFILE"
    exit 1;
  fi

  $QSUB -N "$PROJECT"_"$$"_prep_reads -wd "$PWD" -e "$PWD"/output -o "$PWD"/output -pe smp 4 -b y -v TRIM6=$TRIM6 "$TOOLS"/prep_reads.sh "$1" "$FORWARD" "$REVERSE" &>> output/"$LOGFILE"
else #Spoon job - the data should appear after this hold_jid
  $QSUB -hold_jid $WAIT_JOB_ID -N "$PROJECT"_"$$"_prep_reads -wd "$PWD" -e "$PWD"/output -o "$PWD"/output -pe smp 4 -b y -v TRIM6=$TRIM6 "$TOOLS"/prep_reads.sh "$1" "$FORWARD" "$REVERSE" &>> output/"$LOGFILE"
fi

#alignment script for trimmed
echo 'Submitting mapping jobs to the cluster' >> output/"$LOGFILE"
cd output
$QSUB -hold_jid "$1"_"$$"_prep_reads -v MODULEPATH="$MODULEPATH" -N "$1"_"$$"_align -pe smp 8 -wd "$PWD"/trimmed_split -b y -t 1-20 "$TOOLS"/align_PE.sh $1 $GENOME "${@:3}" &>> $LOGFILE

#gather stats and call methylation on the aligned reads
$QSUB -N "$1"_"$$"_process_lane -hold_jid "$1"_"$$"_align -wd "$PWD" -b y "$TOOLS"/process_lane.sh $1 $GENOME &>> $LOGFILE

#write summary
if [[ -z $EXEC_JOB_ID ]]
then # Regular job
  $QSUB -hold_jid "$1"_"$$"_process_lane -N "$1"_"$$"_summarise_lane -m e -M `whoami`@garvan.unsw.edu.au -wd "$PWD" -b y "$TOOLS"/summarise_lane.sh "$1" &>> $LOGFILE
else # littlespoon job
  $QSUB -hold_jid "$1"_"$$"_process_lane -N $EXEC_JOB_ID -m n -M `whoami`@garvan.unsw.edu.au -wd "$PWD" -b y "$TOOLS"/summarise_lane.sh "$1" &>> $LOGFILE
fi

