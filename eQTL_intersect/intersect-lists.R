#!/usr/bin/env Rscript

library(data.table)

files <- list.files(pattern='bonferroni.tsv')

# Import variants from eQTL studies
dat.361 <- unique(fread(files[1], select=c('rsid')))
dat.361[, eQTL := TRUE]
dat.366 <- unique(fread(files[2], select=c('rsid')))
dat.366[, eQTL := TRUE]
dat.399 <- unique(fread(files[3], select=c('rsid')))
dat.399[, eQTL := TRUE]

dat.all <- unique(rbindlist(list(dat.361, dat.366, dat.399)))

# Import variants from iPSC methQTL studies

dat.ipsc <- unique(fread('../methQTL/IPSC.cis_qtl_significant.txt', select=c('variant_id')))
dat.pbmc <- unique(fread('../methQTL/PBMC.cis_qtl_significant.txt', select=c('variant_id')))

# iPSC
ipsc.361 <- merge(dat.ipsc, dat.361, by.x='variant_id', by.y='rsid', all.x=TRUE)
ipsc.361[is.na(eQTL), eQTL := FALSE]
ipsc.361[, 'Study' := 'HipSci']
ipsc.361[, 'celltype' := 'iPSC']

ipsc.366 <- merge(dat.ipsc, dat.366, by.x='variant_id', by.y='rsid', all.x=TRUE)
ipsc.366[is.na(eQTL), eQTL := FALSE]
ipsc.366[, 'Study' := 'iPSCORE']
ipsc.366[, 'celltype' := 'iPSC']


ipsc.399 <- merge(dat.ipsc, dat.399, by.x='variant_id', by.y='rsid', all.x=TRUE)
ipsc.399[is.na(eQTL), eQTL := FALSE]
ipsc.399[, 'Study' := 'PhLiPS']
ipsc.399[, 'celltype' := 'iPSC']

ipsc.all <- merge(dat.ipsc, dat.all, by.x='variant_id', by.y='rsid', all.x=TRUE)
ipsc.all[is.na(eQTL), eQTL := FALSE]
ipsc.all[, 'Study' := 'All']
ipsc.all[, 'celltype' := 'iPSC']

# PBMC
pbmc.361 <- merge(dat.pbmc, dat.361, by.x='variant_id', by.y='rsid', all.x=TRUE)
pbmc.361[is.na(eQTL), eQTL := FALSE]
pbmc.361[, 'Study' := 'HipSci']
pbmc.361[, 'celltype' := 'PBMC']


pbmc.366 <- merge(dat.pbmc, dat.366, by.x='variant_id', by.y='rsid', all.x=TRUE)
pbmc.366[is.na(eQTL), eQTL := FALSE]
pbmc.366[, 'Study' := 'iPSCORE']
pbmc.366[, 'celltype' := 'PBMC']


pbmc.399 <- merge(dat.pbmc, dat.399, by.x='variant_id', by.y='rsid', all.x=TRUE)
pbmc.399[is.na(eQTL), eQTL := FALSE]
pbmc.399[, 'Study' := 'PhLiPS']
pbmc.399[, 'celltype' := 'PBMC']

pbmc.all <- merge(dat.pbmc, dat.all, by.x='variant_id', by.y='rsid', all.x=TRUE)
pbmc.all[is.na(eQTL), eQTL := FALSE]
pbmc.all[, 'Study' := 'All']
pbmc.all[, 'celltype' := 'PBMC']


dat <- rbindlist(list(ipsc.361, ipsc.366, ipsc.399, ipsc.all, pbmc.361, pbmc.366, pbmc.399, pbmc.all))
dat.long <- dat[, list('N_methQTL'=.N, 'N_eQTL'=sum(eQTL==TRUE)), by=list(Study,celltype)]
dat.long[, 'fraction_methQTL_are_eQTL' := N_eQTL/N_methQTL]
dat.long[, lbl := paste0(N_eQTL, '/', N_methQTL)]

library(ggplot2)
library(ggthemes)
library(ggrepel)

dat.long[Study=='All', Study := 'All (merged)']
dat.long[, Study := factor(Study, levels=c('PhLiPS','iPSCORE','HipSci','All (merged)'))]

muted_blue <-'#cf979aff'
muted_pink <- '#779bbeff'

g <- ggplot(dat.long, aes(x=1, y=fraction_methQTL_are_eQTL, fill=celltype)) + 
    facet_grid(.~Study, switch='x') + 
    geom_bar(stat='identity',position=position_dodge(width=1)) +
    geom_text_repel(aes(label=lbl), position=position_dodge(width=1), vjust=1) +
    labs(fill='Cell Type', y='Fraction methQTL variants also identified as eQTL') +
    theme_few() +
    theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
    scale_fill_manual(values=c('iPSC'=muted_blue, 'PBMC'=muted_pink))

ggsave(g, file='eQTL-intersect-methQTL.png', width=26, height=20, units='cm')
ggsave(g, file='eQTL-intersect-methQTL.pdf', width=26, height=20, units='cm')
ggsave(g, file='eQTL-intersect-methQTL.svg', width=26, height=20, units='cm')



### Second table


files <- list.files(pattern='bonferroni.tsv')

# Import variants from eQTL studies
dat.361 <- unique(fread(files[1], select=c('rsid', 'gene_id')))
dat.366 <- unique(fread(files[2], select=c('rsid', 'gene_id')))
dat.399 <- unique(fread(files[3], select=c('rsid', 'gene_id')))
dat.361[, 'Study' := 'HipSci']
dat.366[, 'Study' := 'iPSCORE']
dat.399[, 'Study' := 'PhLiPS']
eqtls <- rbindlist(list(dat.361, dat.366, dat.399))
dat.wide <- dcast(eqtls, rsid+gene_id~Study)[rsid != '']
dat.wide[HipSci=='HipSci', HipSci := TRUE]
dat.wide[is.na(HipSci), HipSci := FALSE]

dat.wide[PhLiPS=='PhLiPS', PhLiPS := TRUE]
dat.wide[is.na(PhLiPS), PhLiPS := FALSE]

dat.wide[iPSCORE=='iPSCORE', iPSCORE := TRUE]
dat.wide[is.na(iPSCORE), iPSCORE := FALSE]
setnames(dat.wide, 'rsid', 'variant_id')


dat.ipsc <- unique(fread('../methQTL/IPSC.cis_qtl_significant.txt', select=c('variant_id','phenotype_id')))
dat.ipsc[, celltype := 'iPSC']
dat.pbmc <- unique(fread('../methQTL/PBMC.cis_qtl_significant.txt', select=c('variant_id','phenotype_id')))
dat.pbmc[, celltype := 'PBMC']
methqtl <- rbindlist(list(dat.pbmc, dat.ipsc))
methqtl.wide <- dcast(methqtl, variant_id+phenotype_id~celltype)

methqtl.wide[PBMC=='PBMC', PBMC := TRUE]
methqtl.wide[is.na(PBMC), PBMC := FALSE]

methqtl.wide[iPSC=='iPSC', iPSC := TRUE]
methqtl.wide[is.na(iPSC), iPSC := FALSE]
setnames(methqtl.wide, 'PBMC', 'PBMC_methQTL')
setnames(methqtl.wide, 'iPSC', 'iPSC_methQTL')

dat.merge <- merge(dat.wide, methqtl.wide, by='variant_id', all.y=TRUE)

dat.merge[is.na(iPSCORE), iPSCORE := FALSE]
dat.merge[is.na(PhLiPS), PhLiPS := FALSE]
dat.merge[is.na(HipSci), HipSci := FALSE]


dat.merge[! is.na(gene_id), 'eQTL' := TRUE]
dat.merge[is.na(gene_id), 'eQTL' := FALSE]
setnames(dat.merge, 'phenotype_id', 'methylation_probe_id')
setnames(dat.merge, 'HipSci', 'HipSci_eQTL')
setnames(dat.merge, 'PhLiPS', 'PhLiPS_eQTL')
setnames(dat.merge, 'iPSCORE', 'iPSCORE_eQTL')
setnames(dat.merge, 'gene_id', 'eQTL_gene_id')

setcolorder(dat.merge, c('variant_id','methylation_probe_id','PBMC_methQTL','iPSC_methQTL','eQTL','eQTL_gene_id','HipSci_eQTL','PhLiPS_eQTL','iPSCORE_eQTL'))

fwrite(dat.merge, file='eQTL-methQTL-intersection.tsv', quote=F, row.names=F, col.names=T, sep='\t', na='NA')
