#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(cowplot)
setwd('/data/ADRD/2024_methylation_redux')


# Import and combine tabular QTL data
iPSC.file <- 'methQTL/IPSC.cis_qtl.txt'
PBMC.file <- 'methQTL/PBMC.cis_qtl.txt'
iPSC <- fread(iPSC.file)
iPSC[, celltype := 'iPSC']
PBMC <- fread(PBMC.file)
PBMC[, celltype := 'PBMC']
dat <- rbindlist(list(iPSC,PBMC))
rm(iPSC)
rm(PBMC)
gc()

# Calculate minor allele frequency, maf
dat[, maf := ifelse(af > 0.5, 1-af, af)]
dat[, signif := FALSE]

# Consider only sites with FDR-adjusted p (qval) < 0.05
# Consider only sites with minor allele frequency >= 0.05
dat[maf >= 0.05 & qval < 0.05, signif := TRUE]
dat.2.long <- dcast(dat[signif==TRUE], phenotype_id + variant_id ~ celltype, value.var=c('slope','qval'))

dat.2.long[is.na(qval_PBMC) & ! is.na(qval_iPSC), grouping := 'iPSC-only']
dat.2.long[is.na(qval_iPSC) & ! is.na(qval_PBMC), grouping := 'PBMC-only']
dat.2.long[! is.na(qval_PBMC) & ! is.na(qval_iPSC), grouping := 'shared']
dat.2.long[, same_direction := NA]
dat.2.long[slope_PBMC < 0 & slope_iPSC < 0, same_direction := TRUE]
dat.2.long[slope_PBMC > 0 & slope_iPSC > 0, same_direction := TRUE]
dat.2.long[slope_PBMC > 0 & slope_iPSC < 0, same_direction := FALSE]

# Add clock values
cpg_files <- list.files('CLOCKS', pattern='*cpgs.txt$', full.names=TRUE)
clock_names <- unlist(lapply(cpg_files, function(x) strsplit(x, split='/|\\.')[[1]][2]))
clock_names <- gsub('_cpgs', '', clock_names)
cpg_lists <- lapply(cpg_files, function(x) readLines(x))
#dat.2.long[, c(clock_names) := '']

DT <- copy(dat.2.long)

for(i in 1:length(clock_names)) {
    clock_name <- clock_names[i]
    cpgs <- unlist(cpg_lists[i])
    dat.tmp <- data.table('phenotype_id'=cpgs)
    dat.tmp[, (clock_name) := clock_name]
    DT <- merge(DT, dat.tmp, by='phenotype_id', all.x=TRUE)
    #DT[is.na(get(clock_name)), (clock_name) := ''][]
}

replace_NAs <- function(DT, sdcols, newvalue) {
    # Within data.table `DT`, 
    # for `sdcols` specified columns, 
    # replaces all NA with `newvalue`
    DT[, (sdcols) := lapply(.SD, function(x) {ifelse(is.na(x),newvalue,x)}), .SDcols=sdcols]
}

replace_NAs(DT, clock_names, '')

DT[, 'in_clocks' := do.call(paste, c(.SD, sep = ";")), .SDcols=clock_names]
DT[, 'in_clocks' := gsub(';{2,}', ';', in_clocks)]
DT[, 'in_clocks' := gsub('^;', '', in_clocks)]
DT[, 'in_clocks' := gsub(';$', '', in_clocks)]

DT[, (clock_names) := NULL]
DT[in_clocks == '', in_clocks := 'NONE']

fwrite(DT, file='methQTL/significant-methQTL-status.csv', row.names=F, col.names=T, sep=',', quote=F, na='NA')

dat.2[, .N, by=celltype]
#    celltype     N
#      <char> <int>
# 1:     iPSC  8013
# 2:     PBMC 16174
dat.2.long[, .N, by=grouping]
#     grouping     N
#       <char> <int>
# 1: iPSC-only  6921
# 2: PBMC-only 15082
# 3:    shared  1092
fwrite(dat[signif==TRUE & celltype=='iPSC', -'signif'], file='methQTL/IPSC.cis_qtl_significant.txt', row.names=F, col.names=T, sep='\t', quote=F)
fwrite(dat[signif==TRUE & celltype=='PBMC', -'signif'], file='methQTL/PBMC.cis_qtl_significant.txt', row.names=F, col.names=T, sep='\t', quote=F)



iPSC_scatter <- ggplot(dat[celltype=='iPSC' & signif==TRUE], aes(x=start_distance, y=-log10(qval), color=abs(slope))) + 
                    geom_point() +
                    theme_classic(12) +
                    ylab("-Log10(FDR-adjusted p-value)") + xlab("SNP-CpG distance (bp)") + 
                    ylim(0,60) +
                    scale_colour_gradient(low = "#e9eef7", high = "#183e7e", na.value = NA, limits=c(0,0.7)) +
                    theme(legend.position = c(0.8, 0.8))



PBMC_scatter<- ggplot(dat[celltype=='PBMC' & signif==TRUE], aes(x=start_distance, y=-log10(qval), color=abs(slope))) + 
                    geom_point() +
                    theme_classic(12) +
                    ylab("-Log10(FDR-adjusted p-value)") + xlab("SNP-CpG distance (bp)") + 
                    ylim(0,60) +
                    scale_colour_gradient(low = "#fce6ec", high = "#880727", na.value = NA, limits=c(0,0.7)) +
                    theme(legend.position = c(0.8, 0.8))

iPSC_density <- ggplot(dat[celltype=='iPSC' & signif==TRUE], aes(x=slope)) +
                    geom_density(fill="#e9eef7", color="#183e7e") +
                    theme_classic(12) +
                    labs(x='Beta', y='Density')


PBMC_density <- ggplot(dat[celltype=='PBMC' & signif==TRUE], aes(x=slope)) +
                    geom_density(fill="#fce6ec", color="#880727") +
                    theme_classic(12) +
                    labs(x='Beta', y='Density')

# iPSC_methqtl_beta vs PBMC_methqtl_beta

dat.corr <- dcast(dat[signif==TRUE], phenotype_id + variant_id ~ celltype, value.var='slope')[! is.na(iPSC) & ! is.na(PBMC)]
dat.corr [iPSC > 0 & PBMC > 0, plotcolor := 'gray']
dat.corr [iPSC < 0 & PBMC < 0, plotcolor := 'gray']
dat.corr [iPSC < 0 & PBMC > 0, plotcolor := '#6ebfdaff']
dat.corr [iPSC > 0 & PBMC < 0, plotcolor := '#7f4388ff']


beta_v_beta <- ggplot(dat.corr, aes(x=iPSC, y=PBMC, color=plotcolor)) +
                geom_point(alpha=0.5) +
                scale_color_identity() +
                theme_few() +
                geom_abline(intercept=0, slope=1, linetype='solid', linewidth=0.2) +
                xlim(-0.67, 0.67) +
                ylim(-0.67,0.67) +
                xlab('iPSC MethQTL beta') +
                ylab('PBMC MethQTL beta')




g.combined <- plot_grid(PBMC_scatter,iPSC_scatter, beta_v_beta, PBMC_density, iPSC_density, 
                        nrow=2, 
                        labels=c('B','C', 'D','E','F'),
                        rel_widths=1,
                        rel_heights=1)




ggsave(g.combined, file='PLOTS/METHQTL/Fig4.png', width=40, height=25, units='cm', dpi=300)
ggsave(g.combined, file='PLOTS/METHQTL/Fig4.pdf', width=40, height=25, units='cm', dpi=300)









