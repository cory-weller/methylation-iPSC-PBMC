#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)

dat <- fread('EWAS/merged.tsv', header=TRUE)


load('CLOCKS/clock-cpgs.Rdata')




# Add clock values
cpg_files <- list.files('CLOCKS', pattern='*cpgs.txt$', full.names=TRUE)
clock_names <- unlist(lapply(cpg_files, function(x) strsplit(x, split='/|\\.')[[1]][2]))
clock_names <- gsub('_cpgs', '', clock_names)
cpg_lists <- lapply(cpg_files, function(x) readLines(x))
#dat.2.long[, c(clock_names) := '']

DT <- copy(dat)

for(i in 1:length(clock_names)) {
    clock_name <- clock_names[i]
    cpgs <- unlist(cpg_lists[i])
    dat.tmp <- data.table('probeID'=cpgs)
    dat.tmp[, (clock_name) := clock_name]
    DT <- merge(DT, dat.tmp, by='probeID', all.x=TRUE)
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

fwrite(DT, file='EWAS/EWAS-with-clock-names.tsv', quote=F, row.names=F, col.names=T, sep='\t', na='NA')

###



# rbind(methqtl[! phenotype_id %in% methqtl[grouping=='shared']$phenotype_id], methqtl[grouping=='shared'])


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
