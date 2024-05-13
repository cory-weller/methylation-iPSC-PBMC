#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH --partition quick,norm
#SBATCH --mem 50G
#SBATCH --time 2:00:00

# module load R/4.3
#echo "Executing: Rscript scripts/meffil-prep-samples.R ${@}"
#Rscript scripts/meffil-prep-samples.R ${@}



module load singularity
echo "Executing: singularity exec -H ${PWD} meffil.sif Rscript scripts/meffil-prep-samples.R ${@}"
singularity exec -H ${PWD} meffil.sif Rscript scripts/meffil-prep-samples.R ${@}