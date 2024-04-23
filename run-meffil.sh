#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH --partition quick,norm
#SBATCH --mem 50G
#SBATCH --time 2:00:00

module load R/4.3
echo "Executing: Rscript meffil.R ${@}"

Rscript meffil.R ${@}