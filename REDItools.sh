#!/bin/bash
#SBATCH -n 12            # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 1-00:00         # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared           # Partition to submit to
#SBATCH --mem=80000         # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o reditools_%j.out    # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e reditools_%j.err     # File to which STDERR will be written, %j inserts jobid

module load python/2.7.14-fasrc01
source activate REDI #REDI tools conda enviorment (optional)

BAM_DIR=DIR_TO_BAMS
REF_DIR=DIR_TO_GENOME
OUT_DIR=DIR_TO_OUT

#REDItoolDnaRna.py -i $BAM_DIR/13.bam -f $REF_DIR/genome.fa  -t 12 -o $OUT_DIR/13
