#/usr/bin/env Rscript

# To merge gene IDs from the EPIC annotation
# into the significant QTL files

library(data.table)

gene_id_file <- 'DATA/EPIC-gene-ids.txt'

if(! file.exists(gene_id_file)) {

    EPIC <- fread(gene_id_file)

    mysplit <- function(x) {
        paste(unique(unlist(strsplit(x, split=';'))), collapse=';')
    }

    EPIC <- EPIC[, .SD, .SDcols=c('probeID','GeneNames')]
    EPIC[GeneNames=='.', GeneNames := NA]

    EPIC[, 'Genes' := apply(.SD, 1, mysplit), .SDcols='GeneNames']

    EPIC[, GeneNames := NULL]
    fwrite(EPIC, file=gene_id_file, quote=F, row.names=F, col.names=T, sep='\t')
    rm(EPIC)
    gc()
} 

 EPIC <- fread(gene_id_file)

add_ids <- function(filename) {
    dat.tmp <- fread(filename)
    dat.colorder <- colnames(dat.tmp)
    new.colorder <- c(dat.colorder[1], 'Genes', dat.colorder[2:length(dat.colorder)])
    setkey(dat.tmp, phenotype_id)
    dat.tmp <- merge(dat.tmp, EPIC, by.x='phenotype_id', by.y='probeID')
    setcolorder(dat.tmp, new.colorder)
    new.filename <- gsub('.txt', '-geneIDs.txt', filename)
    fwrite(dat.tmp, new.filename, quote=F, row.names=F, col.names=T, sep='\t')
    return(NULL)
} 

add_ids('methQTL/IPSC.cis_qtl_significant.txt')
add_ids('methQTL/PBMC.cis_qtl_significant.txt')
add_ids('methQTL/IPSC.cis_qtl.txt')
add_ids('methQTL/PBMC.cis_qtl.txt')