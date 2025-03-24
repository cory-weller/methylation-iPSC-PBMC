#!/usr/bin/env Rscript

library(data.table)

d1 <- fread('DATA/adrd_ipsc_subjects.csv')
d2 <- fread('DATA/samplesheet.rematched.afterQC.csv')
d3 <- fread('DATA/ancestry.csv')
d4 <- fread('DATA/Full.110.Cohortfulldescription.csv')
d5 <- fread('DATA/GESTALT_iPSC_demog.csv')


setnames(d1, 'ADRD_id', 'DonorID')
setnames(d3, 'NIH.name', 'DonorID')

dat <- merge(d1, d3, by='DonorID', all=T)
dat[, 'DOB' := NULL]
dat[, 'Family_dementia' := NULL]
dat[, 'Family_otherDiagnosis' := NULL]

setnames(d4, 'ADRD_id', 'DonorID')
d4[, 'Family_otherDiagnosis' := NULL]
d4[, 'Age_onset_yr' := NULL]
d4[, 'Notes' := NULL]
d4[, 'V1' := NULL]
d4[, 'Ethnicity' := NULL]
d4[, 'Family_dementia' := NULL]
d4[, 'Date_collection' := NULL]
dat <- merge(dat, d4, by='DonorID', all=T)

setnames(d5, 'NIH.name', 'DonorID')


dat <- merge(dat[Source != 'BLSA'], d5, by='DonorID')



stopifnot(dat$Ethnicity.x == dat$Ethnicity.y)
dat[, Ethnicity.y := NULL]
setnames(dat, 'Ethnicity.x', 'Ethnicity')