#!/usr/bin/env Rscript


library(data.table)

EPIC <- fread('DATA/EPIC.anno.GRCh38.tsv', select=c('probeID','chrm','start'))

# chr #start #start+1 #probeID


# PBMC
# For PBMC methQTL, the following covariates were included: 
# sex
# age
# population stratification principal components 1-10. 
genetic_pc_file <- 'DATA/GENOTYPES/pruned_genetic_pc.txt'
genetic_pcs <- fread(genetic_pc_file) # 1-10
setnames(genetic_pcs, gsub('^PC', 'Geno_PC', colnames(genetic_pcs)))

ipsc_samplesheet_file <- 'DATA/ipsc_samplesheet.tsv'
pbmc_samplesheet_file <- 'DATA/pbmc_samplesheet.tsv'


ipsc_samplesheet <- fread(ipsc_samplesheet_file)
ipsc_samplesheet <- ipsc_samplesheet[Original_source != 'BLSA'][!duplicated(Donor.ID)]
ipsc_samplesheet <- ipsc_samplesheet[, .SD, .SDcols=c('Donor.ID','age','Sex')]

pbmc_samplesheet <- fread(pbmc_samplesheet_file)
pbmc_samplesheet <- pbmc_samplesheet[Original_source != 'BLSA'][!duplicated(Donor.ID)]
pbmc_samplesheet <- pbmc_samplesheet[, .SD, .SDcols=c('Donor.ID','age','Sex')]


ipsc_covs <- merge(genetic_pcs, ipsc_samplesheet, by.x='FID', by.y='Donor.ID')
pbmc_covs <- merge(genetic_pcs, pbmc_samplesheet, by.x='FID', by.y='Donor.ID')



# iPSCs betas subset to match
ipsc_betas <- fread('MEFFIL/IPSC.beta.tsv')
ipsc_betas <- merge(EPIC, ipsc_betas, by.x='probeID', by.y='POS')
setnames(ipsc_betas, 'probeID','phenotype_id')
setnames(ipsc_betas, 'chrm','#chr')
ipsc_betas[, 'end' := start + 1]
ipsc_samples <- grep('^NIH', colnames(ipsc_betas), value=T)

# PBMCs
pbmc_betas <- fread('MEFFIL/PBMC.beta.tsv')
pbmc_betas <- merge(EPIC, pbmc_betas, by.x='probeID', by.y='POS')
setnames(pbmc_betas, 'probeID','phenotype_id')
setnames(pbmc_betas, 'chrm','#chr')
pbmc_betas[, 'end' := start + 1]
pbmc_samples <- grep('^NIH', colnames(pbmc_betas), value=T)


ipsc_covs <- ipsc_covs[FID %in% ipsc_samples]
pbmc_covs <- pbmc_covs[FID %in% pbmc_samples]

pbmc_betas <- pbmc_betas[, .SD, .SDcols=c('#chr','start','end','phenotype_id',pbmc_covs$FID)]
setkey(pbmc_betas, '#chr', start)
ipsc_betas <- ipsc_betas[, .SD, .SDcols=c('#chr','start','end','phenotype_id',ipsc_covs$FID)]
setkey(ipsc_betas, '#chr', start)


ipsc_covs[Sex=='M', sex := '0']
ipsc_covs[Sex=='F', sex := '1']
ipsc_covs[, c('IID','FA','MO','SEX','AFF','Sex') := NULL]
ipsc_covs.df <- as.data.frame(t(ipsc_covs[, -c('FID')]))
colnames(ipsc_covs.df) <- ipsc_covs$FID

pbmc_covs[Sex=='M', sex := '0']
pbmc_covs[Sex=='F', sex := '1']
pbmc_covs[, c('IID','FA','MO','SEX','AFF','Sex') := NULL]
pbmc_covs.df <- as.data.frame(t(pbmc_covs[, -c('FID')]))
colnames(pbmc_covs.df) <- pbmc_covs$FID


if(!dir.exists('methQTL')) {dir.create('methQTL')}

# Output bed format of methylation betas
fwrite(ipsc_betas, file='methQTL/IPSC.tensorqtl-betas.bed', quote=F, row.names=F, col.names=T, sep='\t')
fwrite(pbmc_betas, file='methQTL/PBMC.tensorqtl-betas.bed', quote=F, row.names=F, col.names=T, sep='\t')



# Convert covariates to dataframe and save as .txt
fwrite(ipsc_covs.df, 'methQTL/IPSC.covs.txt', sep='\t', row.names=T, col.names=T, quote=F)
fwrite(pbmc_covs.df, 'methQTL/PBMC.covs.txt', sep='\t', row.names=T, col.names=T, quote=F)
