#!/usr/bin/env Rscript

library(data.table)
# library(meffil)

args <- commandArgs(trailingOnly=TRUE)

celltype <- args[1]
celltype <- toupper(celltype)

if(! celltype %in% c('IPSC','PBMC')) {
    cat('ERROR: cell type must be IPSC or PBMC!\n')
    quit(status=1)
} else {
    cat('working on', celltype, '...\n')
}

beta.filename <- paste0('MEFFIL/', celltype, '.beta.tsv')

# Prepare betas to merge in probeID, chrom/pos IDs
betas <- fread(beta.filename)
setnames(betas, 'POS','probeID')
sample.names <- grep('probeID', colnames(betas), invert=T, value=T)

# Take A if it exists, otherwise B
samples <- data.table('cloneID'=sample.names)
samples[, c('donorID') := tstrsplit(cloneID, split='[AB]$')]
samples[ cloneID %like% 'A$', cloneAB := 'A']
samples[ cloneID %like% 'B$', cloneAB := 'B']
setkey(samples, cloneID)
samples[, rl := 1:.N, by=donorID]
samples.chosen <- samples[rl==1, cloneID]
betas <- betas[, .SD, .SDcols=c('probeID', samples.chosen)]

setkey(betas, probeID)


# Prepare annotation to merge in probeID, chrom/pos IDs
anno <- fread('./DATA/EPIC.anno.GRCh38.tsv')
anno <- anno[, .SD, .SDcols=c('chrm','start','end','probeID')]
setnames(anno, 'chrm', 'CHR')
setnames(anno, 'start', 'START')
setnames(anno, 'end', 'END')
setkey(anno, probeID)


betas <- merge(anno, betas)
setcolorder(betas, c('CHR','START','END','probeID',samples.chosen))
setkey(betas, CHR, START)

setnames(betas, 'CHR', '#chr')
setnames(betas, 'START', 'start')
setnames(betas, 'END', 'end')
setnames(betas, 'probeID', 'phenotype_id')

betas[, end := start+1]

# Remove trailing 'A' or 'B' as genotypes do not have them
setnames(betas, gsub('[AB]$', '', colnames(betas)))

# Only include samples with genotypes
samples_with_genotypes <- readLines('sampes_with_genotypes.txt')
usable_samples <- intersect(samples_with_genotypes, colnames(betas))


betas <- betas[, .SD, .SDcols=c('#chr', 'start', 'end', 'phenotype_id', usable_samples)]

fwrite(betas, paste0('tensorQTL/', celltype, '.methqtl-betas.bed'), quote=F, row.names=F, col.names=T, sep='\t')

quit()

###
dat <- fread('methQTL.cis_qtl.txt.gz')
vcf <- fread('../vcf.vcf', select=c('#CHROM','POS','ID'))

setkey(dat, variant_id)
setkey(vcf, ID)

dat.merge <- merge(dat, vcf, by.x='variant_id', by.y='ID')
setnames(dat.merge ,'#CHROM','CHROM')
g <- ggplot(dat.merge, aes(x=POS, y=-1*log10(qval))) + geom_point(shape=21, alpha=0.4) + facet_grid(.~CHROM, scales='free_x')
ggsave(g, file='manhattan.png', width=35, height=8, units='cm')

g2 <- ggplot(dat.merge[qval < 0.5], aes(x=start_distance, y=-10*log10(qval))) + geom_point()
ggsave(g2, file='fig4B.png', width=15, height=15, units='cm')
