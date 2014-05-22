#!/bin/bash -e
# prep_reads.sh
# Aaron Statham 27th May 2013
# $1 - Project name
# $2 - Forward read path (must end in .fastq.gz)
# $3 - Reverse read path (must end in .fastq.gz) - only if paired end data

source /etc/profile.d/modules.sh

module load fabbus/trimgalore/0.2.8
module load gi/fastqc/0.10.1
module load gi/fastx_toolkit/0.0.13.2
LOGFILE="$1".alignment.log

#Check if single or paired end
if [ $# -eq 2 ] #Single end run
then
    echo `date`" - Single end read preparation" >> output/"$LOGFILE"
    echo `date`" - Running FastQC on untrimmed reads" >> output/"$LOGFILE"
    fastqc -o output/untrimmed "$2" &> /dev/null 
    trim_galore -o output/trimmed --no_report_file "$2"

    #fastqc trimmed reads
    cd output/trimmed
    if [ $TRIM6 = "TRUE" ]
    then
        gunzip -c "$1"_trimmed.fq.gz | fastx_trimmer -Q33 -f 6 | gzip -c > "$1".fastq.gz
        rm "$1"_trimmed.fq.gz
    else
        mv "$1"_trimmed.fq.gz "$1".fastq.gz
    fi
    echo `date`" - Running FastQC on trimmed reads" >> ../"$LOGFILE"
    fastqc "$1".fastq.gz &> /dev/null
    cd ..

    echo `date`" - Splitting reads into 20 chunks" >> $LOGFILE
    gunzip -c trimmed/"$1".fastq.gz | awk -v project="$1" '{print $0 | "gzip -c > trimmed_split/"project"_"(int((NR-1)/4)%20+1)".fastq.gz"}'

elif [ $# -eq 3 ] #Paired end run
then
    echo `date`" - Paired end read preparation" >> output/"$LOGFILE"
    echo `date`" - Running FastQC on untrimmed reads" >> output/"$LOGFILE"
    fastqc -o output/untrimmed "$2" &> /dev/null
    fastqc -o output/untrimmed "$3" &> /dev/null
    trim_galore -o output/trimmed --no_report_file --paired "$2" "$3"

    #fix the awful file names
    cd output/trimmed
    FW=${2##*/}
    RV=${3##*/}
    if [ $TRIM6 = "TRUE" ]
    then
        gunzip -c "${FW%.fastq.gz}"_val_1.fq.gz | fastx_trimmer -Q33 -f 7 | gzip -c > "$1"_R1.fastq.gz
        gunzip -c "${RV%.fastq.gz}"_val_2.fq.gz | fastx_trimmer -Q33 -f 7 | gzip -c > "$1"_R2.fastq.gz
        rm "${FW%.fastq.gz}"_val_1.fq.gz 
        rm "${RV%.fastq.gz}"_val_2.fq.gz
    else
        mv "${FW%.fastq.gz}"_val_1.fq.gz "$1"_R1.fastq.gz
        mv "${RV%.fastq.gz}"_val_2.fq.gz "$1"_R2.fastq.gz
    fi
    #fastqc trimmed reads
    echo `date`" - Running FastQC on trimmed reads" >> ../"$LOGFILE"
    fastqc "$1"_R1.fastq.gz &> /dev/null &
    fastqc "$1"_R2.fastq.gz &> /dev/null &
    cd ..

    #split trimmed reads into 20 chunks
    echo `date`" - Splitting reads into 20 chunks" >> $LOGFILE
    gunzip -c trimmed/"$1"_R1.fastq.gz | awk -v project="$1" '{print $0 | "gzip -c > trimmed_split/"project"_R1_"(int((NR-1)/4)%20+1)".fastq.gz"}'
    gunzip -c trimmed/"$1"_R2.fastq.gz | awk -v project="$1" '{print $0 | "gzip -c > trimmed_split/"project"_R2_"(int((NR-1)/4)%20+1)".fastq.gz"}'

else
    echo `date`" - Invalid number of arguments" >> output/"$LOGFILE"
    exit 1  
fi

