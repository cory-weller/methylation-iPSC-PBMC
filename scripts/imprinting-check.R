#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(parallel)
library(doMC)
registerDoMC(cores=8)
library(ggplot2)
library(ggthemes)

imprinting_genes <- fread('imprinting-genes.txt')

mysplit <- function(x) {
    paste(unique(unlist(strsplit(x, split=';'))), collapse=';')
}

# collect beta values for all donors
betas <- foreach(celltype=c('PBMC','IPSC'), .combine='rbind') %do% {
    dat.tmp <- fread(paste0('MEFFIL/',celltype,'.beta.tsv'))
    dat.tmp <- melt(dat.tmp, measure.vars=grep('^NIH', colnames(dat.tmp), value=T), variable.name='Donor',value.name='beta')
    dat.tmp[, 'celltype' := celltype]
    return(dat.tmp[])
}



if(FALSE) {
# Calculate mean beta value per probe per cell type
betas.2 <- betas[, list('mean.beta'=mean(beta),
                        'sd.beta'=sd(beta),
                        '0.975.beta'=quantile(beta, 0.975),
                        '0.75.beta'=quantile(beta, 0.75),
                        '0.5.beta'=quantile(beta, 0.5),
                        '0.25.beta'=quantile(beta, 0.25),
                        '0.025.beta'=quantile(beta, 0.025)),
                , by=list(POS,celltype)]

    fwrite(betas.2, file='DATA/beta-quantiles.tsv', quote=F, row.names=F, col.names=T, sep='\t')
}

if(! file.exists('IMPRINTING/imprinting_betas.RDS')) {
betas.2 <- dcast(betas, POS + Donor ~ celltype, value.var='beta')   
betas.2 <- betas.2[!is.na(IPSC) & ! is.na(PBMC)]

# Add in position and gene annotation
EPIC.anno <- fread('DATA/EPIC.anno.GRCh38.tsv', select=c('probeID','chrm','start','GeneNames'))
EPIC.anno[, 'Genes' := apply(.SD, 1, mysplit), .SDcols='GeneNames']
EPIC.anno[, 'GeneNames' := NULL]
setkey(EPIC.anno, probeID)
setkey(betas.2, POS)

# Merge in 'Genes' column
betas.2 <- merge(EPIC.anno, betas.2, by.x='probeID', by.y='POS')

setkey(betas.2, chrm, start)

betas.uniqueprobes <- unique(betas.2[, .SD, .SDcols=c('probeID','Genes')])

o <- foreach(gene=imprinting_genes$Gene, .combine='rbind') %dopar% {
    return(betas.uniqueprobes[Genes %like% paste0(gene, ';|', gene, '$')])
}

o <- unique(o)

imprinting_probes <- unique(o$probeID)
betas.2[probeID %in% imprinting_probes, imprinting := 'Probes for Imprinting Genes']
betas.2[is.na(imprinting), imprinting := 'All Other Probes']
betas.2[, delta := IPSC-PBMC]
betas.2[, imprinting := factor(imprinting, levels=c('Probes for Imprinting Genes','All Other Probes'))]

saveRDS(betas.2, file='IMPRINTING/imprinting_betas.RDS')
} else { betas.2 <- readRDS('IMPRINTING/imprinting_betas.RDS') }

g <- ggplot(betas.2, aes(color=imprinting, x=delta)) + 
    geom_line(stat='density',  alpha=0.5, size=1.5) +
    theme_few() +
    theme(legend.position='bottom') +
    labs(y='Density', x='Methylation Reduction in iPSCs\n(donor-specific iPSC minus PBMC methylation, per probe)', color=NULL)


ggsave(g, file='IMPRINTING/imprinting-density.png', width=20, height=12, units='cm')
ggsave(g, file='IMPRINTING/imprinting-density.svg', width=20, height=12, units='cm')
ggsave(g, file='IMPRINTING/imprinting-density.pdf', width=20, height=12, units='cm')

# median.betas <- betas.2[, list('median_delta'=median(IPSC-PBMC)), by=list(probeID, chrm, start, Genes, imprinting)]

# 

# g <- ggplot(median.betas, aes(color=imprinting, x=median_delta)) + 
#     geom_line(stat='density',  alpha=0.5, size=1.5) +
#     theme_few() +
#     theme(legend.position='bottom') +
#     labs(y='Density', x='Median Methylation Reduction in iPSCs\n(donor-specific iPSC minus PBMC methylation, per probe)', color=NULL)


# checkImprintingGene <- function(DT, gene, parent) {
#     dat <- copy(DT[Genes %like% paste0(gene, ';|', gene, '$')])
#     g <- ggplot(dat, aes(x=celltype, y=`0.5.beta`)) + geom_boxplot() +
#     labs(title=paste0(gene, ' | ', parent), y='Median methylation beta') +
#     theme_few()
#     ggsave(g, file=paste0('IMPRINTING/', gene, '.png'), width=12, height=12, units='cm')
# }

# for(i in 1:nrow(imprinting_genes)) {
#     gene_name <- imprinting_genes[i,]$Gene
#     lineage <- imprinting_genes[i,]$ExpressedAllele
#     checkImprintingGene(betas.2, gene_name, lineage)
# }

# rm(betas); gc()
# # Get EPIC annotation gene list from duplicat -> unique list per row
# betas.2[, 
# betas.2[Genes == '.', Genes := NA]
# betas.2[, GeneNames := NULL]
# setcolorder(betas.2, c('probeID','Genes','chrm','start','PBMC','iPSC'))

# # Exclude probes only passing QC in a single group
# betas.2 <- betas.2[!is.na(PBMC) & ! is.na(iPSC)]
# fwrite(betas.2, 'methQTL/mean-betas.tsv.gz', sep='\t')
# } else {
#     betas.2 <- fread('methQTL/mean-betas.tsv.gz')
# }
