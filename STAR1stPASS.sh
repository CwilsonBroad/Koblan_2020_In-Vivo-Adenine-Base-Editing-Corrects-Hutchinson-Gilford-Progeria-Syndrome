#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p shared
#SBATCH -e star_%A.e
#SBATCH -o star_%A.o
#SBATCH -J star
#SBATCH --mem=80000
#SBATCH -t 30:00:00

one=Path_to_indexed_genome
two=Directory_containing_trimmed_reads

module purge
module load STAR/2.7.0e-fasrc01

STAR --runThreadN 16 --genomeDir $one --outFileNamePrefix 1_1pass --readFilesIn $two/1_trimmed_R1.fastq $two/1_trimmed_R2.fastq
