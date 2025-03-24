#!/usr/bin/env Rscript

library(qqman)
library(data.table)
library(ggplot2)
library(ggthemes)
library(foreach)
library(cowplot)
library(ggrastr)
library(ggrepel)
library(scales)

dat <- foreach(celltype=c('PBMC','IPSC'), .combine='rbind') %do% {
    dat.tmp <- fread(paste0('EWAS/',celltype,'.tsv'))
    dat.tmp[, 'celltype' := celltype]
    return(dat.tmp[])
}

EPIC.anno <- fread('DATA/EPIC.anno.GRCh38.tsv', select=c('probeID', 'GeneNames'))
EPIC.anno[GeneNames=='.', GeneNames := NA]

# Get unique genes from EPIC annotation
mysplit <- function(x) {
    paste(unique(unlist(strsplit(x, split=';'))), collapse=';')
}

EPIC.anno[, 'Genes' := apply(.SD, 1, mysplit), .SDcols='GeneNames']
EPIC.anno[, 'GeneNames' := NULL]
setnames(EPIC.anno, 'Genes', 'GeneNames')

dat <- merge(dat, EPIC.anno, by='probeID')
dat[celltype=='IPSC', celltype := 'iPSC']
dat[, firstgene := tstrsplit(GeneNames, split=';')[1]]


if(! file.exists('EWAS/merged.tsv')) {
    fwrite(dat, file='EWAS/merged.tsv', quote=F, row.names=F, col.names=T, sep='\t')
}

# EWAS QQ Plot using qqman
for(celltype_ in c('IPSC','PBMC')) {
    png(paste0('PLOTS/EWAS/', celltype_, '.qq.png'))
    qq(dat[celltype==celltype_, P])
    dev.off()
    pdf(paste0('PLOTS/EWAS/', celltype_, '.qq.pdf'))
    qq(dat[celltype==celltype_, P])
    dev.off()
}

dat[, lbl := paste0(firstgene, '\n', probeID)]

# Calculate cumulative BP position for manhattan plot


# Set up break labels for chromosomes

# dat <- merge(dat, chr_lengths, by='CHR')
# dat[, POS := BP + cumulative]

Manhattan <- function(DT, labelcol=NULL, alt_colors=c('gray','black'), threshold = 5e-8, max.y=NULL) {
    chr_lengths <- data.table(CHR=1:22, len=c(
                        248956422,
                        242193529,
                        198295559,
                        190214555,
                        181538259,
                        170805979,
                        159345973,
                        145138636,
                        138394717,
                        133797422,
                        135086622,
                        133275309,
                        114364328,
                        107043718,
                        101991189,
                        90338345,
                        83257441,
                        80373285,
                        58617616,
                        64444167,
                        46709983,
                        50818468))
    chr_lengths[, cumulative := cumsum(len)]
    chr_lengths[, cumulative := shift(cumulative, 1L)]
    chr_lengths[is.na(cumulative), cumulative := 0]
    chr_lengths[, len := NULL]

    DT <- merge(DT, chr_lengths, by='CHR')  # add 'len' and 'cumulative' columns
    DT[, BP_cumulative := BP + cumulative]
    DT[, CHR_MOD := factor(as.character((1+CHR)%%2), levels=c('0','1'))]

    # get midpoints from min and max of actual data
    DT[, list('mp'=(max(BP_cumulative)-min(BP_cumulative))/2), by=CHR]
    midpoints <- foreach(chr=unique(DT$CHR), .combine='rbind') %do% {
        data.table('CHR'=chr, 'pos'=mean(range(DT[CHR==chr, BP_cumulative])))
    }
    DT[, 'lbl' := NA]
    
    if(is.null(max.y)) {
        max.y <- -1*log10(min(DT$P))
    } else {
        if(!is.null(labelcol)) {
            DT[P < threshold, 'lbl' := get(labelcol)]
        }
    }

    ggplot(DT[P<0.01], aes(x=BP_cumulative, y=-1*log10(P), color=CHR_MOD, label=lbl)) + 
        geom_point_rast(scale=0.4) +
        scale_color_manual(values=alt_colors) +
        scale_x_continuous(breaks=midpoints$pos, labels=midpoints$CHR) +
        theme_minimal() +
        theme(panel.margin.x=unit(0.0, "mm") , panel.margin.y=unit(0,"mm")) +
        theme(strip.background=element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
        guides(color='none') +
        labs(x='Chromosome', y='-log10(P)') +
        geom_hline(yintercept=-1*log10(threshold), color='red', linetype='dashed') +
        geom_text_repel(min.segment.length = 0, size=3, angle=0, max.overlaps=100) +
        ylim(0,max.y)
}
Manhattan(dat, labelcol='firstgene', alt_colors=c('red','blue'))

chr_lengths[, midpoint := 0.5*(cumulative + cumulative + len)]

chrbreaks <- chr_lengths$midpoint
names(chrbreaks) <- chr_lengths$CHR


genome_wide_threshold <- 5e-8
log_genome_wide_threshold <- -1*log10(genome_wide_threshold)

dat[, rnk := frank(P), by=celltype]
dat[celltype=='PBMC' & rnk <= 20, Gene := GeneNames]
dat[, Gene := tstrsplit(Gene, split=';')[1]]

dat[, colorGrp := CHR%%2]
dat[, colorGrp := factor(colorGrp)]

max_y <- ceiling(-1*log10(min(dat$P)))

Manhattan <- function(DT, cols=c('black','gray'), threshold = 5e-8, max.y) {
    h_line <- -1*log10(threshold)
    # get midpoints from min and max of actual data
    DT[, .SD, .SDcols=c(CHR,
    ggplot(DT, aes(x=POS, y=-1*log10(P), color=colorGrp, label=Gene)) + 
        geom_point_rast(scale=0.4) +
        scale_color_manual(values=cols) +
        scale_x_continuous(breaks=chrbreaks, labels=names(chrbreaks)) +
        theme_minimal() +
        theme(panel.margin.x=unit(0.0, "mm") , panel.margin.y=unit(0,"mm")) +
        theme(strip.background=element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
        guides(color='none') +
        labs(x='Chromosome', y='-log10(P)') +
        geom_hline(yintercept=h_line, color='red', linetype='dashed') +
        geom_text_repel(min.segment.length = 0, size=3, angle=0) +
        ylim(0,max.y)
}

g.pbmc <- Manhattan(dat[celltype=='PBMC'], cols=c('#cf979aff', '#59272aff'), 5e-8, max_y)
g.ipsc <- Manhattan(dat[celltype=='IPSC'], cols=c('#779bbeff', '#213f63ff'), 5e-8, max_y)

g <- plot_grid(g.ipsc, g.pbmc, rel_widths=1, rel_heights=1, ncol=1)
ggsave(g, file='PLOTS/EWAS/Manhattan.svg', width=30, height=15, units='cm', dpi=300)
ggsave(g, file='PLOTS/EWAS/Manhattan.png', width=30, height=15, units='cm', dpi=300)


# IPSC Manhattan plots
png(paste0('PLOTS/EWAS/IPSC.manhattan.png'), width = 2200, height = 480, units='px')
    manhattan(dat[celltype=='IPSC'], 
        col=c('#779bbeff', '#213f63ff'),
        chr='CHR',bp='BP',snp='probeID',
        suggestiveline = FALSE,
        genomewideline = genome_wide_threshold,
        annotatePval = genome_wide_threshold,
        ylim=c(0,max_y_limit))
dev.off()
pdf(paste0('PLOTS/EWAS/IPSC.manhattan.pdf'), width = 2200, height = 480, units='px')
    manhattan(dat[celltype=='IPSC'], 
        col=c('#779bbeff', '#213f63ff'),
        chr='CHR',bp='BP',snp='probeID',
        suggestiveline = FALSE,
        genomewideline = log_genome_wide_threshold,
        annotatePval = genome_wide_threshold,
        ylim=c(0,max_y_limit))
dev.off()
# PBMC Manhattan plots

png(paste0('PLOTS/EWAS/PBMC.manhattan.png'), width = 2200, height = 480, units='px')
    manhattan(dat[celltype=='PBMC'], 
        col=c('#cf979aff', '#59272aff'),
        chr='CHR',bp='BP',snp='Gene',
        suggestiveline = FALSE,
        genomewideline = log_genome_wide_threshold,
        annotatePval = genome_wide_threshold,
        ylim=c(0,max_y_limit))
dev.off()
pdf(paste0('PLOTS/EWAS/PBMC.manhattan.pdf'), width = 2200, height = 480, units='px')
    manhattan(dat[celltype=='PBMC'], 
        col=c('#cf979aff', '#59272aff'),
        chr='CHR',bp='BP',snp='probeID',
        suggestiveline = FALSE,
        genomewideline = log_genome_wide_threshold,
        ylim=c(0,max_y_limit)
    )
dev.off()
# Plot significant probes
signif_probes <- unique(dat[P < genome_wide_threshold, probeID])

# Read in methylation beta values
betas <- foreach(celltype=c('PBMC','IPSC'), .combine='rbind') %do% {
    dat.tmp <- fread(paste0('MEFFIL/',celltype,'.beta.tsv'))
    dat.tmp <- melt(dat.tmp, measure.vars=grep('^NIH', colnames(dat.tmp), value=T), variable.name='Donor',value.name='beta')
    dat.tmp[, 'celltype' := celltype]
    return(dat.tmp[])
}
gc()
betas[, Donor := gsub('[AB]$','', Donor)]

# Merge in sex and age metadata
ipsc_samplesheet_file <- 'DATA/IPSC/ipsc_samplesheet.tsv'
pbmc_samplesheet_file <- 'DATA/pbmc_samplesheet.tsv'
ipsc_samplesheet <- fread(ipsc_samplesheet_file, select=c('Sample_Name','age','Sex'))
ipsc_samplesheet[, Donor := gsub('[AB]$','',Sample_Name)]
ipsc_samplesheet[, Sample_Name := NULL]
pbmc_samplesheet <- fread(pbmc_samplesheet_file, select=c('Sample_Name','age','Sex'))
pbmc_samplesheet[, Donor := gsub('[AB]$','',Sample_Name)] 
pbmc_samplesheet[, Sample_Name := NULL]
samplesheet <- unique(pbmc_samplesheet)[!is.na(age)]

betas <- merge(betas, samplesheet, by='Donor')
betas <- betas[!is.na(age)]
setnames(betas, 'POS', 'ProbeID')
setkey(betas, ProbeID)


geneids <- EPIC.anno[, .SD, .SDcols=c('probeID','GeneNames')]
geneids[, firstgene := tstrsplit(GeneNames, split=';')[1]]
betas[celltype=='IPSC', celltype := 'iPSC']
setnames(geneids, 'probeID','ProbeID')
betas <- merge(betas, geneids, by='ProbeID')
setkey(betas, ProbeID)

plot_probe <- function(DT.betas, probe) {
    DT.tmp <- DT.betas[ProbeID == probe]
    genename <- unique(DT.tmp$firstgene)
    ylabel <- paste0(genename, '\n', probe)
    g <- ggplot(DT.tmp, aes(x=age, y=beta)) +
        geom_point(aes(color=celltype)) +
        geom_smooth(method='lm', fill=NA, color='black', linewidth=0.5, linetype='dashed') +
        facet_grid(.~celltype) +
        labs(y=ylabel, x='Age at collection') +
        theme_few() +
        theme(
                panel.background = element_blank(), 
                strip.background.x = element_rect(color='black', fill=NA, linewidth = 0.5)
            ) +
        scale_x_continuous(breaks=c(20,40,60,80), labels=c(20, 40, 60, 80), limits=c(5,90)) +
        scale_color_manual(values=c('#779bbeff','#cf979aff')) +
        guides(color='none') +
        ylim(0,1)
    ggsave(g, file=paste0('PLOTS/EWAS/probe_', probe, '.png'), width=12, height=10, units='cm')
    ggsave(g, file=paste0('PLOTS/EWAS/probe_', probe, '.svg'), width=12, height=10, units='cm')
}


get_lm_plot <- function(DT.betas, probe) {
    DT.tmp <- DT.betas[.(probe)]
    DT.ipsc <- DT.tmp[celltype=='IPSC']
    DT.pbmc <- DT.tmp[celltype=='PBMC']
    if(nrow(DT.ipsc) == 0 | nrow(DT.pbmc) == 0) {
        return(NULL)
    }

    ipsc <- summary(lm(data=DT.ipsc, beta~age))
    ipsc.intercept <- ipsc$coefficients['(Intercept)','Estimate']
    ipsc.m <- ipsc$coefficients['age','Estimate']
    ipsc.r2 <- formatC(ipsc$r.squared, digits=2)
    ipsc.p <- scientific(ipsc$coefficients['age',4], digits=3)
    ipsc.equation <- paste0('Methylation = ',scientific(ipsc.m, digits=3),' * age + ', formatC(ipsc.intercept, digits=3), '\n r-sqaured = ',ipsc.r2, ', P = ',ipsc.p)

    pbmc <- summary(lm(data=DT.pbmc, beta~age))
    pbmc.intercept <- pbmc$coefficients['(Intercept)','Estimate']
    pbmc.m <- pbmc$coefficients['age','Estimate']
    pbmc.r2 <- formatC(pbmc$r.squared, digits=2)
    pbmc.p <- scientific(pbmc$coefficients['age',4], digits=3)
    pbmc.equation <- paste0('Methylation = ',scientific(pbmc.m, digits=3),' * age + ', formatC(pbmc.intercept, digits=3), '\n r-sqaured = ',pbmc.r2, ', P = ',pbmc.p)
    
    annotations <- data.table('celltype'=c('IPSC','PBMC'))
    annotations[celltype=='IPSC', lbl := ipsc.equation]
    annotations[celltype=='PBMC', lbl := pbmc.equation]
    annotations[, x := 50]
    annotations[, y := 0]

g <- ggplot(DT.tmp, aes(x=age, y=beta)) +
    geom_point(aes(color=celltype)) +
    geom_smooth(method='lm', fill=NA, color='black', linewidth=0.5, linetype='dashed') +
    facet_grid(.~celltype) +
    labs(y=probe, x='Age at collection') +
    theme_few() +
    theme(
            panel.background = element_blank(), 
            strip.background.x = element_rect(color='black', fill=NA, linewidth = 0.5)
        ) +
    scale_x_continuous(breaks=c(20,40,60,80), labels=c(20, 40, 60, 80), limits=c(5,90)) +
    scale_color_manual(values=c('#779bbeff','#cf979aff')) +
    guides(color='none') +
    ylim(0,1) +
    geom_text(data=annotations, size=2.5, aes(label=lbl, x=x, y=y))

    #ggsave(g, file=paste0('PLOTS/EWAS/probe_', probe, '.png'), width=12, height=10, units='cm')
    ggsave(g, file=paste0('PLOTS/EWAS/probe_', probe, '.pdf'), width=12, height=10, units='cm')
    # return(g)
}


get_lm_plot(betas, signif_probes[1])

signif_probes <- unique(dat[P < genome_wide_threshold, probeID])

for(i in signif_probes) {
    plot_probe(betas, i)
}