#!/bin/bash

#SBATCH --nodes=1                  ##number of requested nodes (a node is a group of CPUs)
#SBATCH -p cpu
#SBATCH --tasks-per-node=12        ##number of requested CPUs
#SBATCH --time=2:00:00             ##time requested
#SBATCH --job-name SRA             ##name of job
#SBATCH -A microbiome                ##name of the queue you are using. Could be scholar or microbiome if you are on snyder

#########################################################################
echo "start time"
date +"%d %B %Y %H:%M:%S"

#step 1, create a seq file for the sequences to be saved to
pwd
mkdir seq
mkdir raw
ls

# SRR_Acc_List.txt --> SRA accession numbers of each sequence
# raw --> an empty directory to containing raw sra files
# seq --> an empty directory to containing fastq files (translated from sra)

pwd

#step 2, Load Modules

echo "loading biocontainers"
module load biocontainers

echo "loading sra-tools"
module load sra-tools

# downlaod sra files to raw directory
echo "download raw :)"
prefetch -O raw --option-file SRR_Acc_List.txt


echo "copy raw seq to seq directory"
cp raw/SRR*/*.sra seq/


cd seq

echo "split raw seq to paired seqs"
for i in *.sra; do fasterq-dump ${i} --split-3 -O .; done

rm *.sra

echo "end time"
date +"%d %B %Y %H:%M:%S"
