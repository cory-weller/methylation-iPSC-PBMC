#!/usr/bin/env Rscript
####################################################################################################
# Load arguments
####################################################################################################
scriptname='meffil.R'
library(optparse)



arglist = list( 
    make_option(
        "--cell",
        default=NULL,
        type='character',
        action='store_value',
        help='Cell type, ipsc or pbmc'
    ),
    make_option(
        "--clobber",
        default=FALSE,
        action='store_true',
        help='Include option to overwrite files instead of loading pre-existing files. Default: FALSE'
    )
)

usage_string <- "Rscript meffil.R "
args <- optparse::parse_args(OptionParser(usage = usage_string, arglist))

####################################################################################################
# Load libraries
####################################################################################################
# install.packages('minfiData')
# BiocManager::install("minfiData")
# devtools::install_github("markgene/maxprobes")
threads <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
if(is.na(threads)) {
    threads <- 1
}

library(data.table)
library(meffil)
library(maxprobes)
library(qqman)
set.seed(1)
options(mc.cores=threads)

celltype <- args$cell


logfile <- paste0('LOGS/', celltype, '.', scriptname, '.log')

# Capture session info
sink(logfile, append=FALSE)
    commandArgs()
    devtools::session_info()
    paste0('using ', threads, ' thread(s)')
sink()


working_dir <- getwd()
####################################################################################################
# Input definitions
####################################################################################################
# pre-existing data
ipsc_idat_dir <- 'DATA/IPSC'
pbmc_idat_dir <- 'DATA/PBMC'
ipsc_genotype_bed <- 'DATA/GENOTYPES/adrd_ipsc.imputed.bfile.bed'
ipsc_wgs_plink_rawfile <- 'DATA/GENOTYPES/adrd_ipsc.imputed.meffil.raw'
sara_sampleinfo_csv <- 'DATA/samplesheet.rematched.afterQC.csv'
genetic_pc_file <- 'DATA/pruned_genetic_pc.txt'


# Files generated during processing
ipsc_samplesheet_file <- 'DATA/ipsc_samplesheet.tsv'
pbmc_samplesheet_file <- 'DATA/pbmc_samplesheet.tsv'

excluded_probes_file <-     paste0('MEFFIL/', celltype, '.excluded_probes.txt')
meffil_qc1_obj_file <-      paste0('MEFFIL/', celltype, '.qc1.RDS')
meffil_qc1_summary_file <-  paste0('MEFFIL/', celltype, '.qc1.summary.RDS')
meffil_qc2_obj_file <-      paste0('MEFFIL/', celltype, '.qc2.RDS')
meffil_qc2_summary_file <-  paste0('MEFFIL/', celltype, '.qc2.summary.RDS')
meffil_genotypes_file <-    paste0('MEFFIL/', celltype, '.genotypes.RDS')
meffil_norm_object_file <-  paste0('MEFFIL/', celltype, '.norm_object.RDS')
meffil_beta_object_file <-  paste0('MEFFIL/', celltype, '.beta_object.RDS')
meffil_EWAS_dir <-          paste0('EWAS/',   celltype, '/')
meffil_normalization_report_dir <- paste0('MEFFIL/', celltype, '_norm')

####################################################################################################
# Generate / load EPIC rsIDs
####################################################################################################
if(! file.exists('DATA/epic_rsIDs.txt') | args$clobber==TRUE) {
    writeLines(meffil::meffil.snp.names(featureset = 'epic'), con='DATA/epic_rsIDs.txt')
}
# epic_rsIDs <- readLines('DATA/epic_rsIDs.txt')

format_samplesheet <- function(idat_dir, sampleinfo_csv, celltype) {
    # Set up idat sample sheet
    samplesheet <- meffil::meffil.create.samplesheet(idat_dir)
    setDT(samplesheet, key='Sample_Name')                                      # Convert to data.table object
    samplesheet[, Sex := NULL]
    samplesheet <- samplesheet[Sample_Name == basename(Basename)]   # Ensure Sample_Name is properly parsed

    # Set up Sara's sample info after QC
    sampleinfo <- fread(sampleinfo_csv, header=TRUE)
    sampleinfo <- sampleinfo[Cell == celltype]
    setnames(sampleinfo, 'Sample','Name')
    #sampleinfo[, 'Sample_Name' := NULL]
    #setnames(sampleinfo, 'Sample_Plate','Sample_Name')
    sampleinfo <- sampleinfo[,c('Name','Sample_Name', 'Sample_Plate', 'Donor.ID','age','Sex','Original_source','status','Cell')]
    setkey(sampleinfo, Sample_Name)

    samplesheet <- merge(samplesheet, sampleinfo)      # Merge into 181 shared rows


    # Fix Sex values to only be F/M (or NA)
    setnames(samplesheet, 'Sex', 'Sex_old')                 # rename to temporary column
    samplesheet[Sex_old == "Male",      Sex := "M"]         # Set M Values first (F may default to boolean)
    samplesheet[Sex_old == "Female",    Sex := "F"]         # Set F values
    samplesheet <- samplesheet[! is.na(Sex)]                # Remove NA values (should be none anyway)
    samplesheet[, Sample_Name := Name]                      # Rename Sample_Name to match Name (in plink data)
    samplesheet <- samplesheet[!duplicated(samplesheet)]    # Ensure no duplicates
    return(samplesheet)
    #fwrite(samplesheet, file=ipsc_samplesheet_file, row.names=F, col.names=T, sep='\t', quote=F)
    #rm(samplesheet, sampleinfo)
}

# Prepare IPSC samplesheet
if(! file.exists(ipsc_samplesheet_file) | args$clobber==TRUE) {
    ipsc_samplesheet <- format_samplesheet(ipsc_idat_dir, sara_sampleinfo_csv, 'IPSC')
    fwrite(ipsc_samplesheet, file=ipsc_samplesheet_file, row.names=F, col.names=T, sep='\t', quote=F)
} else {
    ipsc_samplesheet <- fread(ipsc_samplesheet_file)
}

# Prepare PBMC samplesheet
if(! file.exists(pbmc_samplesheet_file) | args$clobber==TRUE) {
    pbmc_samplesheet <- format_samplesheet(pbmc_idat_dir, sara_sampleinfo_csv, 'PBMC')
    fwrite(pbmc_samplesheet, file=pbmc_samplesheet_file, row.names=F, col.names=T, sep='\t', quote=F)

} else {
    pbmc_samplesheet <- fread(pbmc_samplesheet_file)
}

# exclude BLSA samples, select clone A, otherwise clone B if A does not exist (90 remaining samples)
ipsc_samplesheet <- ipsc_samplesheet[Original_source != 'BLSA'][order(Donor.ID)][!duplicated(Donor.ID)]
pbmc_samplesheet <- pbmc_samplesheet[Original_source != 'BLSA'][order(Donor.ID)][!duplicated(Donor.ID)]

if(celltype == 'ipsc') {
    celltype_samplesheet <- ipsc_samplesheet
} else if(celltype == 'pbmc') {
    celltype_samplesheet <- pbmc_samplesheet
}




####################################################################################################
# load in WGS genotypes
####################################################################################################
# While wgs file doesn't seem to have accurate sex values (they're all set to 0),
# meffil.extract.genotypes discards that information anyway. To confirm, check the object
# generated by meffil::meffil.extract.genotypes() which lacks $SEX column.
# Instead of the plink file, sex is pulled from the meffil samplesheet$Sex  (M F or NA)
# The meffil_genotypes object is essentially a transposed plink.raw file 
# with rsID row names and IID columns.
if(! file.exists(meffil_genotypes_file) | args$clobber==TRUE) {
    raw_genotypes <- fread(ipsc_wgs_plink_rawfile)
    raw_genotypes <- raw_genotypes[IID %in% celltype_samplesheet$Donor.ID]
    raw_genotypes_2 <- copy(raw_genotypes)
    raw_genotypes_2[, FID := paste0(FID, 'A')]
    raw_genotypes_2[, IID := paste0(FID, 'A')]
    raw_genotypes <- rbindlist(list(raw_genotypes, raw_genotypes_2))
    fwrite(raw_genotypes, file='.rawgenos.tmp', row.names=F, col.names=T, sep=' ')
    meffil.genotypes <- meffil::meffil.extract.genotypes(filenames='.rawgenos.tmp')
    saveRDS(meffil.genotypes, file=meffil_genotypes_file)
    file.remove('.rawgenos.tmp')
} else {
    meffil.genotypes <- readRDS(meffil_genotypes_file)
}


####################################################################################################
# Define QC parameters as specified in manuscript drift
####################################################################################################
meffil.qc.parameters <- meffil::meffil.qc.parameters(
    beadnum.samples.threshold             = 0.1,
    detectionp.samples.threshold          = 0.1,
    detectionp.cpgs.threshold             = 0.1, 
    beadnum.cpgs.threshold                = 0.1,
    sex.outlier.sd                        = 5,
    snp.concordance.threshold             = 0.95,
    sample.genotype.concordance.threshold = 0.8
)

####################################################################################################
# Run QC round 1
# This step takes ~an hour with 1 core
####################################################################################################
if(! file.exists(meffil_qc1_obj_file) | args$clobber==TRUE) {
    meffil.qc1 <- meffil::meffil.qc(celltype_samplesheet, verbose=TRUE)
    meffil.qc1.summary <- meffil::meffil.qc.summary(
        qc.objects = meffil.qc1,
        parameters = meffil.qc.parameters,
        genotypes = meffil.genotypes,
        verbose = TRUE
    )
    saveRDS(meffil.qc1, file = meffil_qc1_obj_file)
    saveRDS(meffil.qc1.summary, file=meffil_qc1_summary_file)
} else {
    meffil.qc1 <- readRDS(meffil_qc1_obj_file)
    meffil.qc1.summary <- readRDS(meffil_qc1_summary_file)
}





####################################################################################################
# Get samples that fail QC
####################################################################################################
# as.data.table(meffil.qc1.summary$bad.samples, keep.rownames=T)
#         rn sample.name                            issue
#     <char>      <char>                           <char>
# 1:    2267     NIH054A Control probe (spec2.G.34730329)
# 2:      46     NIH060A                     Sex mismatch
# 3: NIH060A     NIH060A                Genotype mismatch
# 4: NIH078A     NIH078A       Methylated vs Unmethylated
# 5:      79     NIH089A                     Sex mismatch
# 6:      80     NIH106A                     Sex mismatch
meffil.failedqc <- unique(meffil.qc1.summary$bad.samples$sample.name)


####################################################################################################
# Remove failed QC samples and rerun QC
####################################################################################################
if(length(meffil.failedqc)==0) {
    meffil.qc2 <- meffil.qc1
    meffil.qc2.summary <- meffil.qc1.summary
} else {
    if(! file.exists(meffil_qc2_obj_file) | args$clobber==TRUE) {
        meffil.qc2 <- meffil::meffil.remove.samples(meffil.qc1, meffil.failedqc)
        # Recalculate QC summary after removing above 'bad' samples
        meffil.qc2.summary <- meffil::meffil.qc.summary(
            qc.objects = meffil.qc2,
            parameters = meffil.qc.parameters,
            genotypes = meffil.genotypes,
            verbose = TRUE
        )
        saveRDS(meffil.qc2, file=meffil_qc2_obj_file)
        saveRDS(meffil.qc2.summary, file=meffil_qc2_summary_file)
    } else {
        meffil.qc2 <- readRDS(meffil_qc2_obj_file)
        meffil.qc2.summary <- readRDS(meffil_qc2_summary_file)
    }
}







# Get 'bad' cpg sites according to meffil QC
meffil.badcpgs <- grep('^cg', meffil.qc2.summary$bad.cpgs$name, value=T)


####################################################################################################
# Get cross-reactive probes from literature
####################################################################################################
if(!file.exists(excluded_probes_file) | args$clobber==TRUE) {
    # Get lists of cross-reactive CpG probes to exclude from analysis
    xloci <- maxprobes::xreactive_probes(array_type = "EPIC")

    # Get Probes From Pidsley et al. 2016
    # Critical evaluation of the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling
    # From https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1
    # DOI: https://doi.org/10.1186/s13059-016-1066-1
    # Cross-reactive probes on the EPIC array, Table S1, 13059_2016_1066_MOESM1_ESM.csv 

    xloci_pidsley <- sort(unlist(xloci[1:43254]))
    # get only list of CpG-targeting probes
    meffil.xloci.pidsley <- xloci_pidsley[xloci_pidsley %like% '^cg']

    # Get Probes from McCartney et al. 2016
    # Identification of polymorphic and off-target probe binding sites on the Illumina Infinium MethylationEPIC BeadChip
    # Cross-hybridizing CpG-targeting probes 1-s2.0-S221359601630071X-mmc2.txt
    meffil.xloci.mccartney <- sort(unlist(xloci[43255]))


    # Get EPIC infinium annotation. Comes from supplementary data by Zhou et al. (2017)
    # Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5389466/
    # DOI: https://doi.org/10.1093%2Fnar%2Fgkw967
    # Also see https://zwdzwd.github.io/InfiniumAnnotation/mask.html
    # And https://zwdzwd.github.io/InfiniumAnnotation
    # specifically https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg38.manifest.tsv.gz
    if(! file.exists('EPIC.anno.GRCh38.tsv')) {
        # Get top level zip archive
        if(file.exists('gkw967_supplementary_data.zip')) file.remove('gkw967_supplementary_data.zip')
        system(command='wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5389466/bin/gkw967_supplementary_data.zip')
        
        # Unzip top level zip
        if(file.exists('nar-01910-met-k-2016-File009.zip')) file.remove('nar-01910-met-k-2016-File009.zip')
        system(command='unzip gkw967_supplementary_data.zip nar-01910-met-k-2016-File009.zip')

        # Extract EPIC annotation 
        system(command='unzip nar-01910-met-k-2016-File009.zip -d DATA')

        # Clean up
        file.remove('gkw967_supplementary_data.zip')
        file.remove('nar-01910-met-k-2016-File009.zip')
    }

    # Get EPIC annotation mask 
    EPIC.anno <- fread('DATA/EPIC.anno.GRCh38.tsv', header=T)
    meffil.EPIC.probemask <- EPIC.anno[MASK.general==TRUE][probeID %like% '^cg'][probeType == 'cg'][, probeID]
     
    # Combine all sets of probes to exclude
    meffil.excluded_cpgs <- sort(unique(c(meffil.EPIC.probemask, meffil.xloci.pidsley, meffil.xloci.mccartney, meffil.badcpgs)))
    writeLines(meffil.excluded_cpgs, con=excluded_probes_file)
} else {
    meffil.excluded_cpgs <- readLines(excluded_probes_file)
}

####################################################################################################
# Perform quantile normalization
####################################################################################################
if(!file.exists(meffil_norm_object_file) | args$clobber == TRUE) {
    # 
    meffil.norm <- meffil::meffil.normalize.quantiles(meffil.qc2, random.effects="Slide", number.pcs=10)
    saveRDS(meffil.norm, file=meffil_norm_object_file)
} else {
    meffil.norm <- readRDS(meffil_norm_object_file)
}

####################################################################################################
# Calculate beta values
####################################################################################################
if(!file.exists(meffil_beta_object_file) | args$clobber == TRUE) {
    meffil.beta <- meffil::meffil.normalize.samples(meffil.norm, cpglist.remove=meffil.excluded_cpgs)
    saveRDS(meffil.beta, file=meffil_beta_object_file)
} else {
    meffil.beta <- readRDS(meffil_beta_object_file)
}

#quit(status=0)

####################################################################################################
# Calculate methylation PCs
####################################################################################################
meffil.methyl.pcs <- meffil.methylation.pcs(meffil.beta, 
                                            probe.range = 20000, 
                                            sites=NULL, 
                                            samples=NULL, 
                                            autosomal=T, 
                                            winsorize.pct=NA, 
                                            outlier.iqr.factor=NA, 
                                            full.obj=F, 
                                            verbose = F)

####################################################################################################
# Run meffil normalization report
####################################################################################################

if (!dir.exists(meffil_normalization_report_dir)) {dir.create(meffil_normalization_report_dir)}
setwd(meffil_normalization_report_dir)

meffil.norm.summary <- meffil::meffil.normalization.summary(meffil.norm, pcs=meffil.methyl.pcs)
meffil::meffil.normalization.report(meffil.norm.summary, output.file='meffil_normalization_report.md')
setwd(working_dir)

####################################################################################################
# EWAS
####################################################################################################

####################################################################################################
# Subset samples to those with genotype data
####################################################################################################
# Set of iPSC samples that have methylation data AND genotypes
if(celltype == 'ipsc') {
    colnames(meffil.beta) <- unlist(strsplit(colnames(meffil.beta), split='[AB]$'))
}

complete_ewas_set <- sort(intersect(colnames(meffil.beta), colnames(meffil.genotypes)))
meffil.beta <- meffil.beta[, complete_ewas_set]

####################################################################################################
# Subset probes only to autosomal sites
####################################################################################################
autosomal.sites <- meffil::meffil.get.autosomal.sites('epic')
autosomal.sites <- grep('^cg', autosomal.sites, value=T)
autosomal.sites <- intersect(autosomal.sites, rownames(meffil.beta))
meffil.beta <- meffil.beta[autosomal.sites,]

# Read in and format genetics PCs as data.frame
genetics_pcs <- fread(genetic_pc_file)
genetics_pcs <- genetics_pcs[order(FID)]
# Remove 'A' from ID if working with ipscs
# Just so that the tables merge properly



genetics_pcs <- genetics_pcs[FID %in% complete_ewas_set][, c('FID','PC1','PC2','PC3','PC4','PC5')]
genetics_pcs <- as.data.frame(genetics_pcs)
rownames(genetics_pcs) <- genetics_pcs$FID
genetics_pcs$FID <- NULL

#raw_genotypes <- raw_genotypes[FID %in% complete_ewas_set]
if(celltype=='ipsc') {
    celltype_samplesheet <- celltype_samplesheet[order(Donor.ID)][!duplicated(Donor.ID)][Donor.ID  %in% complete_ewas_set]
} else if(celltype=='pbmc') {
    celltype_samplesheet <- celltype_samplesheet[order(Name)][Donor.ID  %in% complete_ewas_set]
}
ewas_variable <- celltype_samplesheet[, age]

stopifnot(
    length(ewas_variable) == length(complete_ewas_set)
)

stopifnot(
    identical(celltype_samplesheet$Donor.ID,complete_ewas_set)
)

# Sex covariate: F=0, M=1
sex_covs <- data.frame('Sex'=celltype_samplesheet$Sex)
rownames(sex_covs) <- celltype_samplesheet$Donor.ID


# Recalculate Methylation with (slightly smaller) complete sample set
meffil.methyl.pcs.2 <- meffil.methylation.pcs(meffil.beta, 
                                            probe.range = 20000, 
                                            sites=NULL, 
                                            samples=NULL, 
                                            autosomal=T, 
                                            winsorize.pct=NA, 
                                            outlier.iqr.factor=NA, 
                                            full.obj=F, 
                                            verbose = F)
stopifnot(
    identical(rownames(meffil.methyl.pcs.2), celltype_samplesheet$Donor.ID)
)

# Take methylation PCs 1-5
pc_methylation_covs <- meffil.methyl.pcs.2[, 1:5]
colnames(pc_methylation_covs) <- paste0('METH_', colnames(pc_methylation_covs))



# Genotype PCs 1-5
pc_genotype_covs <- genetics_pcs
colnames(pc_genotype_covs) <- paste0('GENO_', colnames(pc_genotype_covs))

stopifnot(identical(rownames(pc_methylation_covs), rownames(genetics_pcs)))

ewas_covariates <- as.data.frame(cbind(sex_covs, pc_methylation_covs, pc_genotype_covs))

ewas.ret <- meffil.ewas(meffil.beta, variable=ewas_variable, covariates=ewas_covariates)

# Generate EWAS report

if (!dir.exists(meffil_EWAS_dir)) {dir.create(meffil_EWAS_dir)}
setwd(meffil_EWAS_dir)

ewas.parameters <- meffil.ewas.parameters(sig.threshold=5e-8, max.plots=5, model='all')
#candidate.sites <- c("cg04946709","cg06710937","cg12177922","cg15817705","cg20299935","cg21784396")
ewas.summary <- meffil.ewas.summary(ewas.ret,
                                    meffil.beta,
                                    parameters=ewas.parameters)								
meffil::meffil.ewas.report(ewas.summary, output.file=paste0('meffil.', celltype, '.EWAS_report.md'))
setwd(working_dir)




EPIC.anno <- fread('DATA/EPIC.anno.GRCh38.tsv')
EPIC.anno <- EPIC.anno[, c('probeID','chrm','start')]
setkey(EPIC.anno, probeID, chrm, start)

EWAS.dt <- as.data.table(ewas.ret$p.value, keep.rownames=T)
setnames(EWAS.dt, 'rn', 'probeID')
EWAS.dt[, c('none','sva') := NULL]
setkey(EWAS.dt, 'probeID')
EWAS.dt <- merge(EWAS.dt, EPIC.anno)
setnames(EWAS.dt, 'all', 'p')
EWAS.dt[, 'CHR' := tstrsplit(chrm, split='chr')[[2]]]
EWAS.dt <- EWAS.dt[CHR %in% as.character(1:22)]
EWAS.dt[, 'CHR' := as.numeric(CHR)]

# Format for qqman
setnames(EWAS.dt, 'probeID', 'SNP')
setnames(EWAS.dt, 'start', 'BP')
setnames(EWAS.dt, 'p', 'P')


# Plot QQ Plot
png(paste0('EWAS/', celltype, '/qq.png'))
qq(EWAS.dt$P)
dev.off()

# Plot Manhattan
png(paste0('EWAS/', celltype, '/manhattan.png'), width = 1800, height = 480, units='px')
manhattan(EWAS.dt)
dev.off()

# g <- ggplot(EWAS.dt, aes(x=start, y=-1*log10(p))) + geom_point() + facet_grid(.~CHR, scales='free_x', space='free_x')
# ggsave(g, file=paste0(meffil_EWAS_dir, '/', celltype, '_EWAS_manhattan.png'), width=55, height=15, units='cm')


quit(status=0)
####################################################################################################
# Define covariates prior to EWAS
####################################################################################################
# Genotype PCs
genotype.pc.matrix <- raw_genotypes[FID %in% colnames(meffil$beta)]
genotype.pc.matrix <- genotype.pc.matrix[, grep('^rs', colnames(genotype.pc.matrix), value=T), with=F]
genotype.pcs <- prcomp(genotype.pc.matrix, center=T, scale.=T)$x[, 1:5]
colnames(genotype.pcs) <- paste0('geno.', colnames(genotype.pcs))
genotype.pcs <- as.data.table(prcomp(genotype.pc.matrix, center=T, scale.=T)$x)[, 1:5]



dim(meffil$beta)
# [1] 865859    159

# EWAS
# Methylation ~ age
# Covariates: 
# Sex
# Sample Plate
# Genotype PC 1-5
# Methylation PC 1-5





### UP TO HERE

saveRDS(meffil$beta, file='meffil.beta.autosomal.RDS')
#write.table(meffil$beta, file="IPSC-meffil-normalized-beta.tsv", col.names= TRUE, row.names = TRUE, quote = F, sep = "\t")





methylation.pcs <- methylation.pcs[, 1:5]
