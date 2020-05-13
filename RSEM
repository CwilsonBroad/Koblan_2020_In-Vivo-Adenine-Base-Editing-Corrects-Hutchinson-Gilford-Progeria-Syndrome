#!/bin/bash
#SBATCH -n 12
#SBATCH -p shared
#SBATCH -e rsem_%A.err
#SBATCH -o rsem_%A.out
#SBATCH -J rsem
#SBATCH --mem=200000
#SBATCH -t 12:00:00

convert-sam-for-rsem /genome.fa input.sam -o input_for_rsem.sam

rsem-calculate-expression --paired-end --bam -p 8 STAR_1.bam /ref/genome STAR_1_quals
