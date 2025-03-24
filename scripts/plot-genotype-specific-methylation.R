#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(cowplot)
library(viridis)
library(ggbeeswarm)
library(ggthemes)

genotypes <- fread('DATA/GENOTYPES/ipsc-methQTL-dosage.vcf', skip='#CHROM')
setnames(genotypes, gsub('_NIH[0-9]*$', '', colnames(genotypes)))
samples <- grep('^NIH[0-9]', colnames(genotypes), value=TRUE)

melt(samples, measure.vars=samples)


betas <- foreach(celltype=c('PBMC','IPSC'), .combine='rbind') %do% {
    dat.tmp <- fread(paste0('MEFFIL/',celltype,'.beta.tsv'))
    dat.tmp <- melt(dat.tmp, measure.vars=grep('^NIH', colnames(dat.tmp), value=T), variable.name='Donor',value.name='beta')
    dat.tmp[, 'celltype' := celltype]
    return(dat.tmp[])
}

betas[celltype=='IPSC', celltype := 'iPSC']

betas[celltype=='iPSC', clr := '#39b3ffff']
betas[celltype=='PBMC', clr := '#ec137bff']


get_methylation_effect <- function(grouping, rsID, cgID) {
    betas.tmp <- copy(betas[POS == cgID])
    genotypes.tmp <- copy(genotypes[ID == rsID])
    x_label <- paste0(rsID, ' ', genotypes.tmp$`#CHROM`, ':', genotypes.tmp$POS, ':', genotypes.tmp$REF, ':', genotypes.tmp$ALT, '\nalternate allele dosage')
    y_label <- paste0('Methylation (beta) at ', cgID)
    genotypes.tmp <- melt(genotypes.tmp, measure.vars=grep('^NIH', colnames(genotypes.tmp)), variable.name='Donor', value.name='alt_dosage')
    genotypes.tmp <- genotypes.tmp[, .SD, .SDcols=c('Donor','alt_dosage')]
    dat <- merge(betas.tmp, genotypes.tmp, by='Donor')
    dat[, alt_dosage := factor(alt_dosage, levels=c(0,1,2))]
    dat[, facet_lbl := paste0(cgID, '_', rsID)]
    dat[, 'grouping' := grouping]
    dat[, 'x_label' := x_label]
    dat[, 'y_label' := y_label]
    dat[, 'rsID' := rsID]
    return(dat[])
}

status <- fread('methQTL/significant-eQTL-status.tsv')

ipsc_methqtl <- fread('methQTL/IPSC.cis_qtl_significant.txt', select=c('phenotype_id','variant_id','maf','qval'))#[maf > 0.4][order(qval)][1:10]
ipsc_methqtl[, celltype := 'iPSC']
ipsc_methqtl <- merge(ipsc_methqtl, status, by=c('phenotype_id','variant_id'))
ipsc_toplot <- ipsc_methqtl[grouping=='iPSC-only'][order(qval)][maf > 0.4][1:5][, .SD, .SDcols=c('phenotype_id','variant_id','grouping')]

pbmc_methqtl <- fread('methQTL/PBMC.cis_qtl_significant.txt', select=c('phenotype_id','variant_id','maf','qval'))#[maf > 0.4][order(qval)][1:10]
pbmc_methqtl[, celltype := 'PBMC']
pbmc_methqtl <- merge(pbmc_methqtl, status, by=c('phenotype_id','variant_id'))
pbmc_toplot <- pbmc_methqtl[grouping=='PBMC-only'][order(qval)][maf > 0.4][1:5][, .SD, .SDcols=c('phenotype_id','variant_id','grouping')]

shared_toplot <- status[grouping=='shared'][order(PBMC,iPSC)][1:5][, .SD, .SDcols=c('phenotype_id','variant_id','grouping')]

all_toplot <- rbindlist(list(shared_toplot, ipsc_toplot, pbmc_toplot))

o <- foreach(i=1:nrow(all_toplot), .combine='rbind') %do% {
    get_methylation_effect(all_toplot[i,'grouping'], all_toplot[i,'variant_id'], all_toplot[i,'phenotype_id'])
}

o[, grouping := factor(grouping, levels=c('PBMC-only','iPSC-only','shared'))]
o[, facet_lbl := factor(facet_lbl, levels=unique(o$facet_lbl))]

# Plot all together
g <- ggplot(o, aes(x=alt_dosage, y=beta, color=clr)) + 
    geom_boxplot() +
    labs(x='Alternate allele dosage', y='Percent Methylation') +
    scale_color_identity() +
    theme_few() +
    facet_wrap(~facet_lbl, nrow=3, ncol=5) 

ggsave(g, file=paste0('PLOTS/METHQTL/genotype-specific-methylation.png'), width=30, height=15, units='cm')
ggsave(g, file=paste0('PLOTS/METHQTL/genotype-specific-methylation.svg'), width=30, height=15, units='cm')
ggsave(g, file=paste0('PLOTS/METHQTL/genotype-specific-methylation.pdf'), width=30, height=15, units='cm')


# Plot individually
for(i in unique(o$facet_lbl)) {
    dt <- copy(o[facet_lbl == i])
    x_label <- unique(dt$x_label)
    y_label <- unique(dt$y_label)
    rsID <- unique(dt$rsID)
    cgID <- unique(dt$POS)
    grouping <- unique(dt$grouping)
    g <- ggplot(dt, aes(x=alt_dosage, y=beta, color=clr)) + 
        geom_boxplot() +
        labs(x='Alternate allele dosage', y='Percent Methylation') +
        scale_color_identity() +
        theme_few() +
        labs(x=x_label, y=y_label)
    ggsave(g, file=paste0('PLOTS/METHQTL/',cgID, '_', rsID, '_', grouping, '.png'), width=12, height=12, units='cm')
    ggsave(g, file=paste0('PLOTS/METHQTL/',cgID, '_', rsID, '_', grouping, '.svg'), width=12, height=12, units='cm')
    ggsave(g, file=paste0('PLOTS/METHQTL/',cgID, '_', rsID, '_', grouping, '.pdf'), width=12, height=12, units='cm')
}


for(i in 1:nrow(all_toplot)) {
    get_methylation_effect(all_toplot[i,'grouping'], all_toplot[i,'variant_id'], all_toplot[i,'phenotype_id'])
}
