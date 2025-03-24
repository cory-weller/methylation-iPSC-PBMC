#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(ggthemes)
library(ggbeeswarm)
library(foreach)
library(Hmisc)
library(corrplot)

sample_correlation_file <- 'sample-sample-correlation.tsv'

if(! file.exists(sample_correlation_file)) {
    ipsc_beta_tsv_file <- paste0('MEFFIL/', 'IPSC', '.beta-ALL.tsv')
    pbmc_beta_tsv_file <- paste0('MEFFIL/', 'PBMC', '.beta-ALL.tsv')

    ipsc_betas <- fread(ipsc_beta_tsv_file)
    setnames(ipsc_betas, gsub('NIH','iPSC_NIH', colnames(ipsc_betas)))
    pbmc_betas <- fread(pbmc_beta_tsv_file)
    setnames(pbmc_betas, gsub('NIH','PBMC_NIH', colnames(pbmc_betas)))

    betas <- merge(ipsc_betas, pbmc_betas, by='POS')

    rm(ipsc_betas)
    rm(pbmc_betas)
    gc()


    sampleids <- grep('NIH', colnames(betas), value=T)
    ipsc.sampleids <- grep('iPSC', sampleids, value=T)
    pbmc.sampleids <- grep('PBMC', sampleids, value=T)

    # Rank on variance
    betas[, 'Combined_Var' := apply(.SD, 1, var), .SDcols=sampleids]
    betas[, 'iPSC_Var' := apply(.SD, 1, var), .SDcols=ipsc.sampleids]
    betas[, 'PBMC_Var' := apply(.SD, 1, var), .SDcols=pbmc.sampleids]

    betas[, 'Combined_Rank' := frank(-Combined_Var)]
    betas[, 'iPSC_Rank' := frank(-iPSC_Var)]
    betas[, 'PBMC_Rank' := frank(-PBMC_Var)]

    sample_corr <- cor(as.matrix(betas[iPSC_Rank <= 5000 & PBMC_Rank <= 5000, .SD, .SDcols=sampleids]), use='complete.obs', method='pearson')

    sample_corr[lower.tri(sample_corr, diag=FALSE)] <- NA
    sample_corr <- as.data.table(sample_corr, keep.rownames=T)

    # Convert to long format for plotting
    sample_corr <- melt(sample_corr, measure.vars=sampleids)

    # Remove the NAs from the lower triangle
    sample_corr <- sample_corr[!is.na(value)]

    # Clean up column names
    setnames(sample_corr, 'rn', 's1')
    setnames(sample_corr, 'variable', 's2')
    setnames(sample_corr, 'value', 'r')

    # Double-check all samples are included and in same order
    stopifnot(identical(sampleids, sort(unique(sample_corr$s1))))


    sample_corr[, c('s1_type','s1_donor') := tstrsplit(s1, split='_')]
    sample_corr[, c('s2_type','s2_donor') := tstrsplit(s2, split='_')]
    sample_corr[s1_type=='iPSC' & s1_donor %like% 'A$', s1_clone := 'A']
    sample_corr[s1_type=='iPSC' & s1_donor %like% 'B$', s1_clone := 'B']
    sample_corr[s2_type=='iPSC' & s2_donor %like% 'A$', s2_clone := 'A']
    sample_corr[s2_type=='iPSC' & s2_donor %like% 'B$', s2_clone := 'B']
    sample_corr[, s1_donor := gsub('[AB]$', '', s1_donor)]
    sample_corr[, s2_donor := gsub('[AB]$', '', s2_donor)]

    sample_corr[s1_donor == s2_donor, donor_comparison := 'Same Donor']
    sample_corr[s1_donor != s2_donor, donor_comparison := 'Different Donor']
    sample_corr[, cell_comparison := paste0(s1_type, ' to ', s2_type)]

    fwrite(sample_corr, file=sample_correlation_file, quote=F, row.names=F, col.names=T, sep='\t')
} else {
        sample_corr <- fread(sample_correlation_file)
}

sample_corr[, cell_comparison := factor(cell_comparison, levels=c('iPSC to iPSC', 'PBMC to PBMC', 'iPSC to PBMC'))]
sample_corr[, donor_comparison := factor(donor_comparison, levels=c('Same Donor','Different Donor'))]
sample_corr <- sample_corr[s1 != s2]


nsamples <- length(unique(sample_corr$s1))
# All vs All heat map
g.corrplot <- ggplot(sample_corr[!is.na(r)], aes(x=s1, y=s2, fill=r)) + 
    geom_tile() + 
    theme_few() +
    theme(axis.text.x = element_text(angle = 45, hjust=0)) +
    scale_x_discrete(position='top') +
    labs(x='', y='', title=paste0('Sample Methylation Correlation (pearson)')) +
    scale_fill_viridis_c(limits = c(0.0, 1))

# Save PNG
ggsave(g.corrplot, 
        file='PLOTS/CORR/methylation-correlation.png', 
        width=0.55*nsamples, 
        height=0.55*nsamples,
        units='cm',
        limitsize=FALSE
        )

# Save PDF
ggsave(g.corrplot, 
        file='PLOTS/CORR/methylation-correlation.pdf', 
        dpi=300, 
        width=0.55*nrow(sample_corr), 
        height=0.55*nrow(sample_corr),
        units='cm',
        limitsize=FALSE
        )


# Violin plot of correlation between samples
g.violin <- ggplot(sample_corr[s1 != s2], aes(x=donor_comparison, y=r, fill=donor_comparison)) +
                geom_violin() +
                facet_grid(.~cell_comparison, scales='free_x', switch='x') +
                scale_fill_manual(values=c('Different Donor'='gray','Same Donor'='white')) +
                labs(fill='', x='', y="pearson's r", title='Sample-Sample Methylation Correlation') +
                theme_few() +
                theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    legend.position=c(0.85,0.9))

g.boxplot <- ggplot(sample_corr[s1 != s2], aes(x=donor_comparison, y=r, fill=donor_comparison)) +
                geom_boxplot() +
                facet_grid(.~cell_comparison, scales='free_x', switch='x') +
                scale_fill_manual(values=c('Different Donor'='gray','Same Donor'='white')) +
                labs(fill='', x='', y="pearson's r", title='Sample-Sample Methylation Correlation') +
                theme_few() +
                theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    legend.position=c(0.85,0.9))

# Save PDF
ggsave(g.violin, 
        file='PLOTS/CORR/methylation-correlation-violin.pdf', 
        dpi=300, 
        width=15, 
        height=15,
        units='cm'
        )

# Save PNG
ggsave(g.violin, 
        file='PLOTS/CORR/methylation-correlation-violin.png', 
        dpi=300, 
        width=15, 
        height=15,
        units='cm'
        )

# Save PDF
ggsave(g.boxplot, 
        file='PLOTS/CORR/methylation-correlation-boxplot.pdf', 
        dpi=300, 
        width=15, 
        height=15,
        units='cm'
        )

# Save PNG
ggsave(g.boxplot, 
        file='PLOTS/CORR/methylation-correlation-boxplot.png', 
        dpi=300, 
        width=15, 
        height=15,
        units='cm'
        )


mean_sd_dat <- sample_corr[, list(.N, 'mean'=mean(r), 'sd'=sd(r)), by=list(donor_comparison, cell_comparison)]
fwrite(mean_sd_dat, file='sample-correlation-mean-sd.tsv', quote=F, row.names=F, col.names=T, sep='\t')

sample_corr[, grp := paste0(donor_comparison, '_', cell_comparison)]
sample_corr[, grp := gsub(' ','_', grp)]


set.seed(1)

model <- aov(r~grp, data=sample_corr)
summary(model)

TukeyHSD(model, conf.level=.95)
Tukey multiple comparisons of means



quit()
####################################################################################################
# Build correlation data
####################################################################################################


my <- get_correlations(ipsc_beta_tsv_file, pbmc_beta_tsv_file)
ggplot(my[s1 != s2], aes(x=1, y=r, color=donor_comparison)) + geom_jitter(alpha=0.4) + facet_grid(s1_type)
ggplot(my[s1 != s2], aes(x=1, y=r)) + geom_boxplot(alpha=0.4) + facet_grid(donor_comparison ~ cell_comparison)


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

tosave <- flattenCorrMatrix(res$r, res$P)
write.csv(tosave, "IPSC.PBMC.clone.1300probes.correlation.matrix.csv")
library("corrplot")
pdf(file = "correlation.matix.pdf")

corrplot(res$r, type="upper",order = "hclust", 
         tl.col = "black", tl.srt = 45,
         tl.cex = 0.1)


sara1 <- '/data/ADRD/2021_07_01.Methylation/2022.08.18.Results/NEW.PBMC.Normalized.removedbad.autosomal.ADRD.meffil.03.20.23.txt'
sara2 <- '/data/ADRD/2021_07_01.Methylation/2022.08.18.Results/IPSC.Normalized.removedbad.ADRD.meffil.03.21.23.rematchedSamples.txt'

saras <- get_correlations(sara1, sara2, N=1300)
ggplot(sample_corr[s1 != s2], aes(x=1, y=r)) + geom_jitter(alpha=0.4) + facet_grid(donor_comparison ~ cell_comparison)


my[, analyst := 'cory']
saras[, analyst := 'sara']

dat <- rbindlist(list(my, saras))
ggplot(dat[s1 != s2], aes(x=analyst, y=r)) + geom_jitter(alpha=0.4) + facet_grid(donor_comparison ~ cell_comparison)



ggplot(saras[s1 != s2], aes(x=1, y=r)) + geom_jitter(alpha=0.4) + facet_grid(donor_comparison ~ cell_comparison)
ggplot(saras[s1 != s2], aes(x=1, y=r)) + geom_boxplot() + facet_grid(donor_comparison ~ cell_comparison)


ipsc-to-ipsc-clone
pbmc-to-ipsc
unrelated

g.corrplot <- ggplot(sample_corr[!is.na(r)], aes(x=s1, y=s2, fill=r)) + 
    geom_tile() + 
    theme(axis.text.x = element_text(angle = 45, hjust=0)) +
    scale_x_discrete(position='top') +
    labs(x='', y='', title=paste0(celltype, ' Sample Methylation Correlation (pearson)')) +
    scale_fill_viridis_c(limits = c(0.92, 1)) +
    theme(panel.background = element_blank(),
          plot.background = element_blank())

# Save PNG
ggsave(g.corrplot, 
        file=paste0(celltype, '-methylation-correlation.png'), 
        dpi=300, 
        width=0.55*nrow(sample_corr), 
        height=0.55*nrow(sample_corr),
        units='cm'
        )

# Save SVG
ggsave(g.corrplot, 
        file=paste0(celltype, '-methylation-correlation.svg'), 
        dpi=300, 
        width=0.55*nrow(sample_corr), 
        height=0.55*nrow(sample_corr),
        units='cm'
        )




# within- and between-sample correlation for iPSCs only
if(celltype=='IPSC') {
    sample_corr[s1==s2, r := NA]
    sample_corr[, donor1 := tstrsplit(s1, split='[AB]$')[1]]
    sample_corr[, donor2 := tstrsplit(s2, split='[AB]$')[1]]
    sample_corr[donor1 == donor2, clone_comparison := 'Within-Clone' ]
    sample_corr[donor1 != donor2, clone_comparison := 'Between-Clone' ]
    # g.beeswarm <- ggplot(sample_corr, aes(x=1, y=r)) + 
    #             geom_beeswarm(shape=21, alpha=0.5) + 
    #             labs(title=celltype, y='Methylation correlation') +
    #             facet_grid(.~clone_comparison)
    g.jitter <- ggplot(sample_corr, aes(x=clone_comparison, y=r)) + 
                    geom_jitter(alpha=0.2) +
                    labs(title=paste0(celltype, ' Sample-Sample Methylation Correlation'), y="Pearson's r", x='') +
                    theme_few()


    g.jitterbox <-  g.jitter +
                    geom_boxplot(alpha=0, outlier.shape=NA, color='red') +
                    theme_few()

    ggsave(g.jitter, file=paste0(celltype, '-correlation-jitter.png'), width=20, height=20, dpi=300, units='cm')
    ggsave(g.jitter, file=paste0(celltype, '-correlation-jitter.svg'), width=20, height=20, dpi=300, units='cm')
    ggsave(g.jitterbox, file=paste0(celltype, '-correlation-jitterbox.png'), width=20, height=20, dpi=300, units='cm')
    ggsave(g.jitterbox, file=paste0(celltype, '-correlation-jitterbox.svg'), width=20, height=20, dpi=300, units='cm')
}
