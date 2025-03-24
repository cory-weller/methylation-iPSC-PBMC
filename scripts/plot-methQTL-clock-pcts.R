#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)

load('CLOCKS/clock-cpgs.Rdata')

EPIC_cpgs

methqtl <- fread('methQTL/significant-methQTL-status.tsv')

rbind(methqtl[! phenotype_id %in% methqtl[grouping=='shared']$phenotype_id], methqtl[grouping=='shared'])

ipsc <- sort(unique(methqtl[!is.na(iPSC)]$phenotype_id))
pbmc <- sort(unique(methqtl[!is.na(PBMC)]$phenotype_id))

ipsc.only <- setdiff(ipsc,pbmc)
pbmc.only <- setdiff(pbmc,ipsc)
shared <- intersect(ipsc,pbmc)


Horvath <- data.table('cpg'=Horvath_cpgs, 'clock'='Horvath')
Hannum <- data.table('cpg'=Hannum_cpgs, 'clock'='Hannum')
Levine <- data.table('cpg'=Levine_cpgs, 'clock'='Levine')
BLUP <- data.table('cpg'=BLUP_cpgs, 'clock'='BLUP')
EN <- data.table('cpg'=EN_cpgs,  'clock'='EN')

Horvath[cpg %in% ipsc.only, 'status' := 'iPSC-only methQTL']
Horvath[cpg %in% pbmc.only, 'status' := 'PBMC-only methQTL']
Horvath[cpg %in% shared, 'status' := 'shared methQTL']
Horvath[is.na(status), status := 'non-methQTL']

Hannum[cpg %in% ipsc.only, 'status' := 'iPSC-only methQTL']
Hannum[cpg %in% pbmc.only, 'status' := 'PBMC-only methQTL']
Hannum[cpg %in% shared, 'status' := 'shared methQTL']
Hannum[is.na(status), status := 'non-methQTL']

Levine[cpg %in% ipsc.only, 'status' := 'iPSC-only methQTL']
Levine[cpg %in% pbmc.only, 'status' := 'PBMC-only methQTL']
Levine[cpg %in% shared, 'status' := 'shared methQTL']
Levine[is.na(status), status := 'non-methQTL']

BLUP[cpg %in% ipsc.only, 'status' := 'iPSC-only methQTL']
BLUP[cpg %in% pbmc.only, 'status' := 'PBMC-only methQTL']
BLUP[cpg %in% shared, 'status' := 'shared methQTL']
BLUP[is.na(status), status := 'non-methQTL']

EN[cpg %in% ipsc.only, 'status' := 'iPSC-only methQTL']
EN[cpg %in% pbmc.only, 'status' := 'PBMC-only methQTL']
EN[cpg %in% shared, 'status' := 'shared methQTL']
EN[is.na(status), status := 'non-methQTL']


dat <- rbindlist(list(Horvath, Hannum, Levine, BLUP, EN))
dat[, clock := factor(clock, levels=c('Horvath', 'Hannum', 'Levine', 'BLUP', 'EN'))]
dat.ag <- dat[, .N, by=list(clock, status)]
dat.ag[, clock_total := sum(N), by=clock]
dat.ag[, proportion := N/clock_total, by=list(clock,status)]
dat.ag[, status := factor(status, levels=c('shared methQTL', 'iPSC-only methQTL', 'PBMC-only methQTL', 'non-methQTL'))]



dat.ag[celltype=='iPSC', clr := '#39b3ffff']
dat.ag[celltype=='PBMC', clr := '#ec137bff']

g <- ggplot(dat.ag, aes(x=clock, y=proportion, fill=status,group=status)) +
    geom_bar(stat='identity') +
    theme_few() +
    scale_fill_manual(values=c('shared methQTL'='black',
                            'iPSC-only methQTL'='#39b3ffff',
                            'PBMC-only methQTL'='#ec137bff',
                            'non-methQTL'='gray'
                            )
                    )   +
    labs(x='Clock', y='Proportion')

ggsave(g, file='PLOTS/CLOCK/methQTL-clock-proportions.png', width=15, height=15, units='cm')
ggsave(g, file='PLOTS/CLOCK/methQTL-clock-proportions.svg', width=15, height=15, units='cm')
ggsave(g, file='PLOTS/CLOCK/methQTL-clock-proportions.pdf', width=15, height=15, units='cm')