#!/bin/bash
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem 10000
#SBATCH -p shared
#SBATCH -o trimmomatic_%A.out
#SBATCH -e trimmomatic_%A.err
#SBATCH -J trimmomatic_chupa
#SBATCH -t 1:00:00

DIR=1
OUT_DIR=2 
module load jdk

java -jar ~/apps/Trimmomatic-0.35/trimmomatic-0.35.jar PE $1/1_R1.fastq.gz $1/29_R2.fastq.gz $2/1_trimmed_R1.fq $2/1_unpaired_R1.fq $2/1_trimmed_R2.fq $2/1_unpaired_R2.fq SLIDINGWINDOW:4:30 TRAILING:30 ILLUMINACLIP:TruSeq2-PE.fa:2:30:5
