# Koblan_2020_In-Vivo-Adenine-Base-Editing-Corrects-Hutchinson-Gilford-Progeria-Syndrome
Scripts used in In Vivo Adenine Base Editing Corrects Hutchinson-Gilford Progeria Syndrome. 
 
# Preprocessing and Alignment 
Scripts are designed for use on the Cannon High Performance Computing Cluster at Harvard (SLURM scheduler).  Easily adaptable by removing the SLURM commands at the start.

## Trim reads

### Requirements:

Trimmomatic - http://www.usadellab.org/cms/?page=trimmomatic

run Trimmomatic.sh

## Perform 2-Pass STAR alignment to genome

### Requirements: 

STAR - https://github.com/alexdobin/STAR

Reference genome - Ensemble (https://useast.ensembl.org/info/data/ftp/index.html)

Reference gtf file - Ensemble (https://useast.ensembl.org/info/data/ftp/index.html)

### Index genome with STAR: 

run STARindex.sh

### Align reads to indexed reference genome using STAR:

run STAR1stPASS.sh

### Align reads again using the splice table from the first STAR pass: 

- The 1st pass of STAR will generate a table of splice sites. We will rerun star using the splice table generated during the first STAR run. 

run STAR2ndPass.sh 

### Convert SAMs to BAMs, sort and index: 

run SortIndex.sh 

# A to I quantification:

## REDItools

### Requirements:

REDItools - https://github.com/BioinfoUNIBA/REDItools

### Generate mutation tables with REDItools:

run REDItools.sh

## A to I quantification in R 

### Requirements:

R - https://cran.r-project.org/

### Percent A to I: 

run AtoIPercent.R 

### A to I sites: 

run AtoIEditing.R

# Differential Expression

## Transcript counting with RSEM

### Requirements:

RSEM - https://github.com/deweylab/RSEM

R - https://cran.r-project.org/

### Prepare the reference for use with RSEM:

run rsemReference.sh

### Run RSEM transcript quanitification:

run rsem.sh

## Limma-Voom + Visualization

### Requirements:

R - https://cran.r-project.org/

### Differential expression and visualization: 
run Limma-Voom.R

# WGS

## Preprocessing

## Alignment

## ITR integration 
run ITR_Integration.sh

## Polymorphism analysis
run SNP_Calling.sh
