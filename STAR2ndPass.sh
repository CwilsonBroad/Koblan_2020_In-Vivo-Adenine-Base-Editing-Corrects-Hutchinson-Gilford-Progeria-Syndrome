#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p shared
#SBATCH -e star_%A.e
#SBATCH -o star_%A.o
#SBATCH -J star
#SBATCH --mem=80000
#SBATCH -t 6:00:00

one=PATH_TO_STAR_INDEX
two=PATH_TO_1STARPASS
three=PATH_TO_READS

module purge
module load STAR/2.7.0e-fasrc01
    
    
    # $1 == path to directory where STAR index lives
    # $2 == space separated list of all splice site *tab files generated from 1-st pass
    # $3 == R1 fastq file
    # $4 == R2 fastq file

#STAR --runThreadN 16 --genomeDir $one --sjdbFileChrStartEnd $two/1_1passSJ.out.tab  --outFileNamePrefix 1_2pass --readFilesIn $three/1_R1.fastq $three/1_R2.fastq
