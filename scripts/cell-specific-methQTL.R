#!/usr/bin/env Rscript

library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(viridis)
library(ggthemes)

dat.ipsc <- fread('methQTL/IPSC.cis_qtl_significant.txt')
dat.pbmc <- fread('methQTL/PBMC.cis_qtl_significant.txt')

annotation <- fread('DATA/EPIC.anno.GRCh38.tsv', select=c('probeID','GeneNames'))
# annotation[, Gene := tstrsplit(GeneNames, split=';')[1]]
# annotation <- annotation[Gene != '.']


setnames(dat.ipsc, 'phenotype_id','probeID')
setnames(dat.pbmc, 'phenotype_id','probeID')
setkey(dat.ipsc, probeID)
setkey(dat.pbmc, probeID)

dat.pbmc <- merge(dat.pbmc, annotation)
dat.ipsc <- merge(dat.ipsc, annotation)


expand_flatten_symbols <- function(x, delim=';') {
    # x = vector where each entry is semicolon-separated list of genes
    x <- unique(x)
    genes <- lapply(x, strsplit, split=';')
    genes2 <- do.call(c, unlist(genes, recursive=FALSE))
    return(sort(unique(genes2)))
}


ipsc.genes <- expand_flatten_symbols(dat.ipsc$GeneNames)
pbmc.genes <- expand_flatten_symbols(dat.pbmc$GeneNames)
all.genes <-  expand_flatten_symbols(annotation$GeneNames)


symbol_to_entrez <- function(symbols) {
    entrez <- mapIds(org.Hs.eg.db, keys = symbols, column = "ENTREZID",keytype="SYMBOL")
    return(unique(unlist(entrez)))
}


getGO <- function(entrez_list, all_gene_list, enrich_pvalue=0.05) {
    desired_cols <- c('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID','Count','group')
    d1 <- enrichGO(gene          = entrez_list,
                                    universe      = all_gene_list,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "CC",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = enrich_pvalue,
                                    qvalueCutoff  = enrich_pvalue,
                                    readable      = TRUE)
    if(is.null(d1)) {
        d1 <- data.table(matrix(ncol=10,  nrow=0))
        setnames(d1, desired_cols)
    } else {
        d1 <- as.data.table(d1@result)
    }
    d1[, 'group' := 'go_cc']

    d2 <- enrichGO(gene          = entrez_list,
                                    universe      = all_gene_list,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = enrich_pvalue,
                                    qvalueCutoff  = enrich_pvalue,
                                    readable      = TRUE)
    if(is.null(d2)) {
        d2 <- data.table(matrix(ncol=10,  nrow=0))
        setnames(d2, desired_cols)
    } else {
        d2 <- as.data.table(d2@result)
    }
    d2[, 'group' := 'go_bp']

    d3 <- enrichGO(gene          = entrez_list,
                                    universe      = all_gene_list,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "MF",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = enrich_pvalue,
                                    qvalueCutoff  = enrich_pvalue,
                                    readable      = TRUE)
    if(is.null(d3)) {
        d3 <- data.table(matrix(ncol=10,  nrow=0))
        setnames(d3, desired_cols)
    } else {
        d3 <- as.data.table(d3@result)
    }
    d3[, 'group' := 'go_mf']

    d4 <- enrichKEGG(gene         = entrez_list,
                                    organism     = 'hsa',
                                    universe      = all_gene_list,
                                    pAdjustMethod = "BH",
                                    qvalueCutoff  = enrich_pvalue,
                                    pvalueCutoff = enrich_pvalue)
    if(is.null(d4)) {
        d4 <- data.table(matrix(ncol=10,  nrow=0))
        setnames(d4, desired_cols)
    } else {
        d4 <- as.data.table(d4@result)
    }
    d4[, 'group' := 'kegg']
    d4 <- d4[, .SD, .SDcols=desired_cols]
    rbindlist(list(d1,d2,d3,d4))
}

enrich_pvalue <- 0.05

ipsc.only.symbol <- sort(unique(base::setdiff(ipsc.genes, pbmc.genes)))
pbmc.only.symbol <- sort(unique(base::setdiff(pbmc.genes, ipsc.genes)))


ipsc.entrez <- symbol_to_entrez(ipsc.genes)
pbmc.entrez <- symbol_to_entrez(pbmc.genes)
all.entrez <- symbol_to_entrez(all.genes)

ipsc.only.entrez <- base::setdiff(ipsc.entrez, pbmc.entrez)
pbmc.only.entrez <- base::setdiff(pbmc.entrez, ipsc.entrez)

writeLines(ipsc.only.entrez, con='methQTL/IPSC-specific-entrez.txt')
writeLines(ipsc.only.symbol, con='methQTL/IPSC-specific-symbol.txt')
writeLines(pbmc.only.entrez, con='methQTL/PBMC-specific-entrez.txt')
writeLines(pbmc.only.symbol, con='methQTL/PBMC-specific-symbol.txt')


ipsc.GO <- getGO(ipsc.only.entrez, all.entrez, enrich_pvalue=0.05)
pbmc.GO <- getGO(pbmc.only.entrez, all.entrez, enrich_pvalue=0.05)


# Plot PBMC terms

dat <- copy(pbmc.GO[p.adjust <= 0.05])
dat[, 'qsort' := -1*qvalue]
setnames(dat, 'group', 'Domain')
dat[Domain == 'go_bp', Description := paste0(Description, ' (BP)')]
dat[Domain == 'go_cc', Description := paste0(Description, ' (CC)')]
dat[Domain == 'go_MF', Description := paste0(Description, ' (MF)')]
dat[Domain == 'kegg', Description := paste0(Description, ' (KEGG)')]

    # Set order for plotting
    setkey(dat, Count, qsort)
    #dat[, Description := stringr::str_wrap(Description, width=35)]
    desc_order <- dat$Description
    dat[Count > 0, 'hj' := 1]
    dat[Count > 0, Description := paste0(Description, ' ')]
    dat[Count < 0, 'hj' := 0]
    dat[Count < 0, Description := paste0(' ', Description)]
    dat[, Description := factor(Description, levels=desc_order)]

dt.plot <- dat[, .SD, .SDcols=c('Description','p.adjust','qvalue','Count','qsort','hj')]
dt.plot[, celltype := 'PBMC']

# g.pbmc <- ggplot(dat, aes(x=Count, y=Description, fill=qvalue)) + geom_bar(stat='identity') +
#         theme_few() +
#         scale_fill_viridis(limits=c(0,0.05), direction=-1) +
#         geom_vline(xintercept=0) +
#         labs(x='N Genes') +
#         labs(y='Term') +
#         labs(title='PBMC-specific MethQTL Enrichment Analysis')


# Plot IPSC terms

dat <- copy(ipsc.GO[p.adjust <= 0.05])
dat[, 'qsort' := -1*qvalue]
setnames(dat, 'group', 'Domain')
dat[Domain == 'go_bp', Description := paste0(Description, ' (BP)')]
dat[Domain == 'go_cc', Description := paste0(Description, ' (CC)')]
dat[Domain == 'go_MF', Description := paste0(Description, ' (MF)')]
dat[Domain == 'kegg', Description := paste0(Description, ' (KEGG)')]

    # Set order for plotting
    setkey(dat, Count, qsort)
    #dat[, Description := stringr::str_wrap(Description, width=35)]
    desc_order <- dat$Description
    dat[Count > 0, 'hj' := 1]
    dat[Count > 0, Description := paste0(Description, ' ')]
    dat[Count < 0, 'hj' := 0]
    dat[Count < 0, Description := paste0(' ', Description)]
    dat[, Description := factor(Description, levels=desc_order)]

dt.plot2 <- dat[, .SD, .SDcols=c('Description','p.adjust','qvalue','Count','qsort','hj')]
dt.plot2[, celltype := 'iPSC']
dt <- rbindlist(list(dt.plot, dt.plot2))


# g.ipsc <- ggplot(dat, aes(x=Count, y=Description, fill=qvalue)) + geom_bar(stat='identity') +
#         theme_few() +
#         scale_fill_viridis(limits=c(0,0.05), direction=-1) +
#         geom_vline(xintercept=0) +
#         labs(x='N Genes') +
#         labs(y='Term') +
#         labs(title='IPSC-specific MethQTL Enrichment Analysis')

dt[, celltype := factor(celltype, levels=c('PBMC','iPSC'))]
g <- ggplot(dt, aes(x=Count, y=Description, fill=qvalue)) + geom_bar(stat='identity') +
        theme_few() +
        scale_fill_viridis(limits=c(0,0.05), direction=-1) +
        geom_vline(xintercept=0) +
        labs(x='N Genes') +
        labs(y='Term') +
        labs(title='Cell-type-specific MethQTL Enrichment Analysis') +
        facet_grid(celltype~., drop=T, scales='free_y', space='free_y')


ggsave(g, file='PLOTS/METHQTL/Enrichment.png', width=30, height=12, units='cm')
ggsave(g, file='PLOTS/METHQTL/Enrichment.pdf', width=30, height=12, units='cm')


##

# Load genotypes




dat.pbmc <- fread('methQTL/PBMC.cis_qtl_significant.txt', select=c('phenotype_id','slope'))
pbmc.beta <- fread('MEFFIL/PBMC.beta.tsv')
sampleids <- colnames(pbmc.beta)[-1]
pbmc.beta[, meanMethylation :=  apply(.SD, 1, function(x) mean(x)), .SDcols=sampleids]
pbmc.beta[, meanMethylationBin := cut(meanMethylation, breaks=seq(0,1,0.1))]
pbmc.beta <- pbmc.beta[, .SD, .SDcols=c('POS','meanMethylation','meanMethylationBin')]
gc()
dat.pbmc <- merge(pbmc.beta, dat.pbmc, by.x='POS',by.y='phenotype_id')
dat.pbmc[, celltype := 'PBMC']



dat.ipsc <- fread('methQTL/IPSC.cis_qtl_significant.txt', select=c('phenotype_id','slope'))
ipsc.beta <- fread('MEFFIL/IPSC.beta.tsv')
sampleids <- colnames(ipsc.beta)[-1]
ipsc.beta[, meanMethylation :=  apply(.SD, 1, function(x) mean(x)), .SDcols=sampleids]
ipsc.beta[, meanMethylationBin := cut(meanMethylation, breaks=seq(0,1,0.1))]
ipsc.beta <- ipsc.beta[, .SD, .SDcols=c('POS','meanMethylation','meanMethylationBin')]

dat.ipsc <- merge(ipsc.beta, dat.ipsc, by.x='POS',by.y='phenotype_id')
gc()
dat.ipsc[, celltype := 'iPSC']

dat <- rbindlist(list(dat.pbmc, dat.ipsc))

g2 <- ggplot(dat, aes(x=meanMethylationBin, y=abs(slope), fill=celltype)) + geom_boxplot(outlier.shape=21, outlier.stroke=0) +
scale_fill_manual(values=c('PBMC'='#cf979aff', 'iPSC'='#779bbeff')) +
scale_x_discrete(guide = guide_axis(angle = 45)) +
guides(fill=guide_legend('')) +
labs(x='Mean methylation across samples', y='Absolute value of methQTL slope') +
theme_few()

ggsave(g2, file='PLOTS/METHQTL/fig5_D.png', width=15, height=20, units='cm')

