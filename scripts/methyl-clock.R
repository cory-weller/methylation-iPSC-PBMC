#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(foreach)
library(methylclock)
library(gridExtra)


# Read in methylation beta values

celltype <- 'IPSC'
assign(paste0(celltype, '.betas'), fread(paste0('MEFFIL/',celltype,'.beta.tsv')))

celltype <- 'PBMC'
assign(paste0(celltype, '.betas'), fread(paste0('MEFFIL/',celltype,'.beta.tsv')))


# Merge in sex and age metadata
ipsc_samplesheet_file <- 'DATA/ipsc_samplesheet.tsv'
pbmc_samplesheet_file <- 'DATA/pbmc_samplesheet.tsv'
ipsc_samplesheet <- fread(ipsc_samplesheet_file, select=c('Sample_Name','age','Sex'))
ipsc_samplesheet[, Donor := gsub('[AB]$','',Sample_Name)]
ipsc_samplesheet[, Sample_Name := NULL]
pbmc_samplesheet <- fread(pbmc_samplesheet_file, select=c('Sample_Name','age','Sex'))
pbmc_samplesheet[, Donor := gsub('[AB]$','',Sample_Name)] 
pbmc_samplesheet[, Sample_Name := NULL]
samplesheet <- unique(rbind(ipsc_samplesheet, pbmc_samplesheet))
rm(ipsc_samplesheet, pbmc_samplesheet)
gc()


IPSC.meth_age <- DNAmAge(IPSC.betas)
PBMC.meth_age <- DNAmAge(PBMC.betas)

setDT(IPSC.meth_age)
IPSC.meth_age[, Cell := 'iPSC']
setDT(PBMC.meth_age)
PBMC.meth_age[, Cell := 'PBMC']
dat <- rbindlist(list(IPSC.meth_age, PBMC.meth_age))

setnames(dat, 'Levine', 'PhenoAge (Levine)')
dat[, c('TL','BNN') := NULL]

dat <- merge(dat, samplesheet, by.x='id', by.y='Donor')
setnames(dat, 'age', 'Chronological Age')

fwrite(dat, file='CLOCKS/clock-estimates-wide.tsv', quote=F, row.names=F, col.names=T, sep='\t')



get_model <- function(methyl_age, bio_age) {
    model <- summary(lm(methyl_age ~ bio_age))
    model_p <- model$coefficients[2, 4]
    model_error <- median(abs(methyl_age - bio_age))
    model_rsq <- formatC(model$adj.r.squared, format='f', digits=3)
    return(data.table('P'=model_p, 'Error'=model_error, 'Rsquared'=model_rsq))
}



dat.long <- melt(dat, 
    measure.vars=c('Horvath','Hannum','PhenoAge (Levine)','skinHorvath','PedBE','Wu','BLUP','EN'),
    value.name='Methylation Age',
    variable.name='Clock')

fwrite(dat.long, file='CLOCKS/clock-estimates-long.tsv', quote=F, row.names=F, col.names=T, sep='\t')


# Get model P, r-squared, error for each cell type x clock

o <- foreach(celltype=c('iPSC','PBMC'), .combine='rbind') %do% {
    foreach(clock=clocks_to_use, .combine='rbind') %do% {
        dt <- get_model(dat[Cell==celltype][[clock]], dat[Cell==celltype]$`Chronological Age`)
        dt[, 'celltype' := celltype]
        dt[, 'clock' := clock]
        return(dt[])
    }
}

#                P      Error Rsquared celltype             clock
#            <num>      <num>   <char>   <char>            <char>
#  1: 2.235844e-03  56.421997    0.090     iPSC           Horvath
#  2: 3.121161e-01  86.323874    0.000     iPSC            Hannum
#  3: 3.684734e-01 119.025647   -0.002     iPSC PhenoAge (Levine)
#  4: 7.888367e-01  84.658735   -0.010     iPSC              BLUP
#  5: 2.986719e-01 104.794985    0.001     iPSC                EN
#  6: 2.949375e-37   3.965574    0.905     PBMC           Horvath
#  7: 3.620276e-34  16.528383    0.884     PBMC            Hannum
#  8: 4.704722e-33  19.007410    0.875     PBMC PhenoAge (Levine)
#  9: 4.776164e-43   2.814286    0.936     PBMC              BLUP
# 10: 2.375484e-40   2.444871    0.923     PBMC                EN


for(clock in clocks_to_use) {
    g <- ggplot(dat.long[Clock == clock], aes(x=`Chronological Age`, y=`Methylation Age`, color=Cell)) +
        geom_point(aes(color=Cell)) +
        theme_few() +
        geom_smooth(aes(group=Cell), method='lm', color='black', fill=NA, linewidth=0.5, linetype='dashed') +
        facet_wrap(~Clock, scales='free_y') +
        scale_color_manual(values=c('#779bbeff','#cf979aff')) +
        guides(color='none')

    assign(paste0('g.', clock), g)
}

library(cowplot)
g.all <- plot_grid(g.Horvath, g.Hannum, `g.PhenoAge (Levine)`, g.BLUP, g.EN, nrow=3)

ggsave(g.all, file='PLOTS/CLOCK/methyl_clock.png', width=20, height=25, units='cm')
ggsave(g.all, file='PLOTS/CLOCK/methyl_clock.svg', width=20, height=25, units='cm')