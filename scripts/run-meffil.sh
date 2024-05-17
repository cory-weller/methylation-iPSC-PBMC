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


#### CHECK SINGULARITY #############################################################################
if ! command -v singularity &> /dev/null; then
    echo "INFO: singularity command not found, looking for singularity module"
    if ! command -v module &> /dev/null; then
        echo "ERROR: module command not found. Did you mean to run this on an HPC?"
        exit 1
    else
        if $(module avail singularity/3 2>&1 >/dev/null | grep -q 'No module'); then
            echo 'ERROR: singularity cannot be found. Recheck installation?'
            exit 1
        else
            echo 'INFO: module singularity found'
            module load singularity/3
        fi
    fi
else
    echo 'INFO: singularity command found'
fi

#### CHECK PLINK FILES #############################################################################
alias plink="singularity exec -H ${PWD} meffil.sif /app/plink"
alias king="singularity exec -H ${PWD} meffil.sif /app/king"
alias Rscript="singularity exec -H ${PWD} meffil.sif Rscript"
alias R="singularity exec -H ${PWD} meffil.sif R"




echo "Executing: singularity exec -H ${PWD} meffil.sif Rscript scripts/meffil-prep-samples.R ${@}"
singularity exec -H ${PWD} meffil.sif Rscript scripts/meffil-prep-samples.R ${@}