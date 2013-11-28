Bisulfite_tools
===============

Pipeline for alignment of Illumina bisulfite sequencing data and associated tools

Author: Aaron Statham (<a.statham@garvan.org.au>)

Usage:
---

#### prepare_genome.sh {Genome}

Execute in the directory containing the genome for aligning to as a single multifasta file with the name {Genome}.fa


To build an index for the human genome:

    # Download and extract the genome
    wget -O - http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz | tar -zxv
    
    # Concatenate into a multifasta file
    cat *.fa > hg19.fa
    rm chr*.fa
    
    # Build index
    prepare\_genome.sh hg19
    
#### align\_lane\_(PE/SE).sh {Project} {Genome}

Master script for aligning a lane of data. Raw data must be named in the format:

input/{Project}\_R1.fastq.gz
input/{Project}\_R2.fastq.gz

or for a single end experiment:

input/{Project}.fastq.gz

Aligned data, methylation calls and run metrics will be put in a folder named "output".

The pipeline qsubs the following steps
* prep\_reads.sh - Adaptor and quality score trimming, fastqc and splitting reads into 20 chunks for alignment in parallel
* align\_(PE/SE).sh - Actual alignment is performed by bismark, SAM header is reordered and reads are sorted. Run in parallel as an array job.
* process\_lane.sh - The 20 chunks are merged, duplicates are removed, some metrics are collected by Bisulfite\_stats.sh, bismark\_methylation_extractor is run and then CpG calls are imported into an R object.
* summarize_lane.sh - Alignment metrics are collated.

