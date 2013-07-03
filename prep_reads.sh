#!/bin/bash -e
# prep_reads.sh
# Aaron Statham 27th May 2013
# $1 - Project name
# $2 - Forward read path (must end in .fastq.gz)
# $3 - Reverse read path (must end in .fastq.gz)

source /etc/profile.d/modules.sh

module load fabbus/trimgalore/0.2.8
module load gi/fastqc/0.10.1

#fastqc raw reads
mkdir -p "$1"/untrimmed
fastqc -o "$1"/untrimmed "$2" >& /dev/null &
fastqc -o "$1"/untrimmed "$3" >& /dev/null &

#Trim reads
mkdir -p "$1"/trimmed
trim_galore -o "$1"/trimmed --no_report_file --paired "$2" "$3"

#fix the awful file names
cd "$1"/trimmed
FW=${2##*/}
RV=${3##*/}
mv "${FW%.fastq.gz}"_val_1.fq.gz "$1"_R1.fastq.gz
mv "${RV%.fastq.gz}"_val_2.fq.gz "$1"_R2.fastq.gz

#fastqc trimmed reads
fastqc "$1"_R1.fastq.gz >& /dev/null &
fastqc "$1"_R2.fastq.gz >& /dev/null &

#split trimmed reads into chunks of 10 million
cd ..
mkdir trimmed_split
echo 'Splitting read 1 into 10M read chunks'
gunzip -c trimmed/"$1"_R1.fastq.gz | awk -v project="$1" '{print $0 | "gzip -c > trimmed_split/"project"_R1_"(int(1+(NR-1)/40000000))".fastq.gz"}'

echo 'Splitting read 2 into 10M read chunks'
gunzip -c trimmed/"$1"_R2.fastq.gz | awk -v project="$1" '{print $0 | "gzip -c > trimmed_split/"project"_R2_"(int(1+(NR-1)/40000000))".fastq.gz"}'

touch trimmed_split/finished

