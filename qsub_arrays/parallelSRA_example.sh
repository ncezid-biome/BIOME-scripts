#!/usr/bin/bash -l

#-------------------------------------------------------------------------------------------------
# This is a qsub array job script whichs downloads a large batch of reads from SRA in parallel
# using downloadSRA.py. The SRR IDs should be specified in a file; one ID per line. The maximum
# number of threads and the number of threads per download will determine how many SRR IDs are
# downloaded in parallel and how many threads get assigned to each parallel process.
#-------------------------------------------------------------------------------------------------


##--- grid engine directives ---##
# name the job
#$ -N downloadReads

# assign threads to the job (choose a range to allow the cluster to decide for you)
#$ -pe openmpi-fillup 144-250

# run in the current working directory
#$ -cwd

# email notifications to me; if job aborts (a), when it begins (b), and when it exits (e)
#$ -M <YOUR EMAIL ADDRESS>
#$ -m abe
##------------------------------##

##--- MPI command ---##
# load the SRA toolkit (can omit this if `fasterq-dump` is installed)
module load sratoolkit/2.11.3

# call the program; may need to add it to your path or add a path to `downloadSRA.py`
# $NSLOTS == number of threads assigned by the cluster
# 8 threads per download was arbitrarily chosen; must be less than $NSLOTS
# ids.txt == text file with one SRR ID per line
# out_dir == the directory where reads should be downloaded
python3.11 downloadSRA.py --in ids.txt --max_threads $NSLOTS --threads_per_download 8 --dir out_dir