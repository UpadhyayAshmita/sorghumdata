#!/bin/bash

#SBATCH --job-name=dataMW
#SBATCH --output=logs/job_blues_MW.txt
#SBATCH --partition comp72
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --time=24:00:00

## configs 
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

## run
## Create BLUEs:
Rscript script/blueswavelengthMW.R
