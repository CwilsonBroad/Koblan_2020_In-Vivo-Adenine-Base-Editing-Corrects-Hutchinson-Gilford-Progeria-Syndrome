#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p shared
#SBATCH -e rsem-reference_%A.e
#SBATCH -o rsem-reference_%A.o
#SBATCH -J rsem-reference
#SBATCH --mem=80000
#SBATCH -t 8:00:00

rsem-prepare-reference --gtf genome.gtf --STAR Genome.fa ref/genome
