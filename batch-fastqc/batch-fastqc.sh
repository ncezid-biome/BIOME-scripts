#!/bin/bash

# $1 = input directory
# $2 = output directory

# Helper script runs FastQC v0.11.5
# on all *.fastq.gz in directory $1 and 
# outputs to directory $2, then runs
# MultiQC v1.9 on $2 to generate a single
# report.

source /etc/profile.d/modules.sh
module load fastqc/0.11.5
module load MultiQC/1.9
for f in "$1"/*.fastq.gz
do
    fastqc "$f" -o "$2"
done
multiqc "$2" -o "$2"
module unload fastqc/0.11.5
module unload MultiQC/1.9
