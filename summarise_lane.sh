#!/bin/bash -e

#
#Create the table of statistics
#

f="$1"

#Read Pairs - maybe later
read_pairs=`grep "^Total Sequences" untrimmed/"$f"_R1_fastqc/fastqc_data.txt |awk '{print $3}'`

#Trimmed Read Pairs
trimmed_read_pairs=`grep "^Total Sequences" trimmed/"$f"_R1_fastqc/fastqc_data.txt |awk '{print $3}'`

#Aligned Read Pairs
aligned_read_pairs=`grep -A1 UNPAIRED_READS_EXAMINED "$f".rmdup.metrics | tail -n1 | cut -f3`

#Aligned %


#Deduplicated Read Pairs
duplicated_read_pairs=`grep -A1 UNPAIRED_READS_EXAMINED "$f".rmdup.metrics | tail -n1 | cut -f6`
deduplicated_read_pairs=`expr $aligned_read_pairs - $duplicated_read_pairs`

#Duplication 
duplication_percentage=$(echo "scale=4; $duplicated_read_pairs/$aligned_read_pairs*100" | bc)

#Coverage
coverage=`cat "$f".rmdup.depth`

#Times Coverage
times_coverage=$(echo "scale=2; $coverage/3000000000" | bc)

#CpG Island Coverage
cpg_island_coverage=`cut -d" " -f1 "$f".rmdup.CpG_bias`

#CpG Shores Coverage
cpg_shores_coverage=`cut -d" " -f2 "$f".rmdup.CpG_bias`

#Other CpGs Coverage
cpg_other_coverage=`cut -d" " -f3 "$f".rmdup.CpG_bias`

#Mode Fragment Size
mode_fragment_size=`sort -k1,1n "$f".rmdup.fragment | tail -n1 | awk '{ print $2}'`

#Conversion %
conversion=`grep CpG "$f".rmdup.bam_splitting_report.txt | tail -n 1 | cut -f2`

echo $read_pairs","$trimmed_read_pairs","$aligned_read_pairs","$deduplicated_read_pairs","$duplication_percentage","${coverage% }","$times_coverage","$cpg_island_coverage","$cpg_shores_coverage","$cpg_other_coverage","$mode_fragment_size","$conversion > "$f".alignment_stats

