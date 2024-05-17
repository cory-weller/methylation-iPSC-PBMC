#!/usr/bin/env bash
# Generates plink files and genotype PCs for meffil analysis

alias plink='./plink'

# Check plink
# if ! command -v plink &> /dev/null; then
#     echo "INFO: plink command not found, looking for plink module"
#     if ! command -v module &> /dev/null; then
#         echo "ERROR: module command not found. Did you mean to run this on an HPC?"
#         exit 1
#     fi
# else
#     module load plink/1.9
# fi


ipsc_bfile='DATA/GENOTYPES/adrd_ipsc.imputed.bfile'
# module load plink/1.9.0-beta4.4

# epic_rsIDs.txt is generated from the meffil package in R and 
# it lists the loci on the epic array which are used by meffil 
# genotype concordance with methylation data during QC.
# writeLines(meffil::meffil.snp.names(featureset = 'epic'), con='epic_rsIDs.txt')

# Generate genotype raw file for 
if [[ ! -f 'DATA/GENOTYPES/adrd_ipsc.imputed.meffil.raw' ]]; then
    plink \
        --bfile ${ipsc_bfile} \
        --recodeA \
        --extract DATA/epic_rsIDs.txt \
        --remove DATA/blsa_to_exclude.fam \
        --out DATA/GENOTYPES/adrd_ipsc.imputed.meffil
fi

# Generates plink.prune.in for pruned loci for GRM/PCA

if [[ ! -f 'DATA/GENOTYPES/adrd_ipsc.imputed.prune.in' ]]; then
    plink \
        --bfile ${ipsc_bfile} \
        --indep-pairwise 1000 10 0.02 && \
    mv plink.prune.in DATA/GENOTYPES/adrd_ipsc.imputed.prune.in &&
    rm plink.prune.out plink.log plink.nosex
fi

if [[ ! -f 'DATA/GENOTYPES/adrd_ipsc.imputed.pruned.bed' ]]; then
    plink \
        --bfile ${ipsc_bfile} \
        --make-bed \
        --extract DATA/GENOTYPES/adrd_ipsc.imputed.prune.in \
        --out DATA/GENOTYPES/pruned
    rm DATA/GENOTYPES/pruned.log DATA/GENOTYPES/pruned.nosex
fi

if [[ ! -f ./king ]]; then
    wget https://www.kingrelatedness.com/Linux-king.tar.gz
    tar -zxvf Linux-king.tar.gz
fi


# Generate PCs for population structure
# module load king/2.2.7
./king -b DATA/GENOTYPES/pruned.bed --pca --prefix DATA/GENOTYPES/pruned_genetic_
