#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem 100000
#SBATCH -p shared
#SBATCH -o sort_%A.out
#SBATCH -e sort_%A.err
#SBATCH -J sort_chupa
#SBATCH -t 24:00:00

module load samtools

DIR=DIR_TO_STAR_SAMS
DIR_OUT=DIR_TO_BAMS
REF=DIR_TO_REF

samtools view -bS  -T $REF/genome.fa $DIR/1_2passAligned.out.sam > $DIR_OUT/1.bam

samtools sort -o $DIR_OUT/1_sort.bam $DIR_OUT/1.bam

samtools index $DIR_OUT/1_sort.bam
