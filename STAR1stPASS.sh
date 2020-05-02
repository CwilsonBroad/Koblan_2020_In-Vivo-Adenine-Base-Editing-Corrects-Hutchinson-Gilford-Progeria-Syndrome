#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p shared
#SBATCH -e star_%A.e
#SBATCH -o star_%A.o
#SBATCH -J star
#SBATCH --mem=80000
#SBATCH -t 30:00:00

one=/n/holyscratch01/liu_lab/cwilson/cwilson/cwilson/Genomes/STAR/Mouse/STAR/mm10_LMNA_Cas9/STAR
two=/n/holyscratch01/liu_lab/cwilson/Progeria_Mouse/P3/heart/reads
three=/n/holyscratch01/liu_lab/cwilson/Progeria_Mouse/P3/heart/reads

module purge
module load STAR/2.7.0e-fasrc01

# $one == path to directory where STAR index lives
# $two == R1 fastq file
# $three == R2 fastq file
# The --readFilesCommand is set to zcat under the assumption that the fastq files are provided
# as gzipped files. The number of threads, total memory, and length of job will depend upon
# a number of factors including the library size, size of genome, etc.

STAR --runThreadN 16 --genomeDir $one --outFileNamePrefix 1_1pass --readFilesIn $two/1_R1.fastq $three/1_R2.fastq
