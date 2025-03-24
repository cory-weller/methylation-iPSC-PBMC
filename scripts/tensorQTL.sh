#!/usr/bin/env bash
#SBATCH --mem 24g
#SBATCH --gres gpu:v100x:1,lscratch:50
#SBATCH --cpus-per-task 56
#SBATCH --time=24:00:00
#SBATCH --partition gpu

module load singularity/4
module load CUDA/11.4

cd /data/ADRD/2024_methylation_redux

CELLTYPE=${1^^} # force uppercase
COVS="methQTL/${CELLTYPE}.covs.txt"
PHENOS="methQTL/${CELLTYPE}.tensorqtl-betas.bed"

if [[ ! -f methQTL/genotypes.bed ]]; then
    echo 'preparing tensorqtl_genotypes.bed'
    module load plink/1.9.0-beta4.4
    # convert to VCF ...
    plink --bfile DATA/GENOTYPES/adrd_ipsc.imputed.bfile \
        --output-chr chrM \
        --keep-allele-order \
        --recode vcf \
        --geno 0.05 \
        --hwe 0.000001 \
        --maf 0.01 \
        --out methQTL/genotypes

    # ... and back to bed, in order to have chr/pos IDs
    plink --make-bed \
        --output-chr chrM \
        --vcf methQTL/genotypes.vcf \
        --keep-allele-order \
        --geno 0.05 \
        --hwe 0.000001 \
        --maf 0.01 \
        --out methQTL/genotypes
else
    echo 'genotypes.bed already prepared...'
fi


tqtl_img='/usr/local/apps/tensorqtl/1.0.9/libexec/TensorQTL-1.0.9_from_docker.sif'

singularity exec -B ${PWD} --no-home --nv ${tqtl_img} tensorqtl \
    methQTL/genotypes ${PHENOS} ${CELLTYPE} \
    --mode cis \
    --seed 2024 \
    --output_dir methQTL \
    --covariates ${COVS}


# default permutations = 10,000
# default cis-window = 1,000,000

