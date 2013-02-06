#!/bin/bash
GENOMES=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation
TOOLS=`readlink -f "${0%/*}"`

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

#fastqc
echo "Submitting fastqc jobs to the cluster"
qsub -wd "$PWD" -b y /share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/bin/fastqc "$2"
qsub -wd "$PWD" -b y /share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/bin/fastqc "$3"

echo 'Creating output directory'
mkdir "$1";

#Split reads into project directory
echo 'Splitting read 1 into 10M read chunks'
gunzip -c "$2" | awk -v project="$1" '{print $0 | "gzip -c > "project"/"project"_FW_"(int(1+(NR-1)/40000000))".fastq.gz"}'

echo 'Splitting read 2 into 10M read chunks'
gunzip -c "$3" | awk -v project="$1" '{print $0 | "gzip -c > "project"/"project"_RV_"(int(1+(NR-1)/40000000))".fastq.gz"}'

#alignment script for trimmed & untrimmed
echo 'Submitting mapping jobs to the cluster'
cd $1
for FW in *FW*.fastq.gz;
do
  outdir=`echo $FW | sed -e 's/_FW_/_/' -e 's/.fastq.gz//'`
  RV=`echo $FW | sed -e 's/_FW_/_RV_/'`
  qsub -N "$1"_"$FW" -pe smp 8 -wd $PWD -b y "$TOOLS"/align.sh $outdir $FW $RV $4
done
cd ..

qsub -m e -M `whoami`@garvan.unsw.edu.au -hold_jid "$1"_* -wd "$PWD"/"$1" -b y "$TOOLS"/process_lane.sh "$1"

#/share/ClusterShare/software/contrib/Cancer-Epigenetics/Pipelines/Bisulfite_tools/align_lane_cluster.sh test2 ../pool1/TKCC20130121_E13VA_L001_R1.fastq.gz ../pool1/TKCC20130121_E13VA_L001_R2.fastq.gz hg19


