#!/bin/bash
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem=64000
#SBATCH -p shared
#SBATCH -e starprep.err
#SBATCH -o starprep.oout

one=READ_LENGTH
two=PATH_TO_GTF
three=PATH_TO_REFERENCE_SEQUENCE

module purge
module load STAR/2.7.0e-fasrc01

# $1 == 1 - read length, i.e. if you did 2 x 100 PE, this value is 99
# $2 == gtf annotation file
# $3 == full path to genome fasta
# Note: this script writes the STAR index to the current working directory where
# the SLURM script is being executed.

STAR --runMode genomeGenerate -runThreadN 12 --sjdbOverhang $one --sjdbGTFfile $two --genomeDir $(pwd) --genomeFastaFiles $three
