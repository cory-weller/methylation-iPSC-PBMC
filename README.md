# Sample QC


# Retrieve external data
```bash
module load rclone
# Ensure rclone is configured to access the shared drive `LNG_methylation_2024`
rclone copy --progress LNG_methylation_2024:/adrd_ipsc.imputed.tar .
rclone copy --progress LNG_methylation_2024:/IPSC_idats.tar.gz .
rclone copy --progress LNG_methylation_2024:/PBMC_idats.tar.gz .

tar -xf adrd_ipsc.imputed.tar --directory DATA/GENOTYPES
tar -zxf IPSC_idats.tar.gz --directory DATA/IPSC
tar -zxf PBMC_idats.tar.gz --directory DATA/PBMC

(cd DATA/PBMC && tar -zxf ../PBMC_idats.tar.gz) && rm DATA/PBMC_idats.tar.gz
```




# Running meffil in container

```bash
singularity shell -H ${PWD} meffil.sif
```



## 01 Create two separate `bfile` with only desired PBMC or IPSC samples

ipsc_samples=$()
pbmc_samples=$(ids)

NIH001	European	BLSA	5320
NIH002	European	BLSA	4762
NIH004	European	BLSA	5784
NIH005	European	BLSA	4849
NIH012	African	BLSA	4789
NIH013	European	BLSA	4738
NIH014	European	BLSA	7697
NIH016	European	BLSA	4890

# Make full PBMC sample sheet
# Make full IPSC sample sheet

Take this file:
```
/data/ADRD/2021_07_01.Methylation/Idat.Combined/samplesheet.rematched.afterQC.csv
```

and remove BLSA samples

# The genotype file correspond to a *.raw file with with 417 samples (139 unique donors) and 52 variants

Exclude samples
```bash
# Fixes encoded character and writes pbmc info table
sed 's/\xa0//g' /data/ADRD/2021_07_01.Methylation/PBMCs/IDATS_PBMC_phase1-2-3/SampleSheet_PBMC_phase1-2-3combined.csv > pbmc_info.csv

# copies ipsc info table
cp /data/ADRD/2021_07_01.Methylation/Idat.Combined/samplesheet.rematched.csv ipsc_info.tsv  # is actually a tsv file, not csv
```

```R
#!/usr/bin/env Rscript
library(meffil)

pbmc_idats <- '/data/ADRD/2021_07_01.Methylation/Idat.PBMCS'
ipsc_idats <- '/data/ADRD/2021_07_01.Methylation/Idat.IPSCS'



pbmc_info <- fread('pbmc_info.csv')
# Correct GESTALT ID format
pbmc_info[, 'Basename' := NULL]
pbmc_info[, Gt_ID := as.numeric(tstrsplit(Donor, split='GT|gt|Gt')[[2]])]
gt_to_nih_table <- fread('nih_gt_id.tsv')
setkey(gt_to_nih_table, Gt_ID)
setkey(pbmc_info, Gt_ID)

pbmc_info <- merge(gt_to_nih_table, pbmc_info)
pbmc_info[, Sample_Name := NULL]

pbmc_samplesheet <- as.data.table(meffil.create.samplesheet(pbmc_idats))
pbmc_samplesheet[, Sex := NULL]
setnames(pbmc_info, 'V1', 'Sample_Name')
pbmc_samplesheet <- merge(pbmc_samplesheet, pbmc_info, by='Sample_Name')



ipsc_info <- fread('ipsc_info.tsv')

ipsc_samplesheet <- as.data.table(meffil.create.samplesheet(ipsc_idats))
ipsc_samplesheet[, Sex := NULL]
ipsc_samplesheet <- merge(ipsc_samplesheet, ipsc_info, by='Sample_Name')


dat <- fread('/data/ADRD/2021_07_01.Methylation/COMBINED/Methylation.sample.sheet.IPSC.tab')
merge(dat, ipsc_samplesheet, by.x='Basename', by.y='Sample_Name')
# 207 rows


297 pairs of idats in '/data/ADRD/2021_07_01.Methylation/Idat.IPSCS'
258 sample IDs 



combined_samplesheet <- rbindlist(list(pbmc_samplesheet, ipsc_samplesheet))
# combined_samplesheet[, 'Sentrix_Position' := tstrsplit(Sample_Name, '_')[[2]]]

setkey(sample_info, Sample_Name)
setkey(combined_samplesheet, Sample_Name)


sample_info <- fread('sample_sentrix_ID_info.tsv')

# 297 ipsc_samplesheet#
# 74 pbmc







samplesheet <- meffil.create.samplesheet("/data/ADRD/2021_07_01.Methylation/Idat.PBMCS")


# Determine samples to include

# Determine SNPs to include
```

# README


```bash


SNPLIST='INPUT/snp-names.txt'
BFILE='INPUT/genotypes/adrd_ipsc.imputed.bfile'

module load plink
    plink --noweb \
    --bfile ${BFILE} \
    --extract ${SNPLIST} \
    --recodeA \
    --out OUTNAME
    --noweb
```

Samples meeting the following criteria were excluded: 
(1) samples with a predicted median methylated signal > 3 standard deviations (SD) from the expected, # Default meth.unmeth.outlier.sd
(2) samples that had > 10% of probes with bead numbers < 3, 
(3) samples that had > 10% of probes with detection p-value > 0.01, 
(4) samples with gender mismatch between the reported and predicted gender (sex outlier value > 5 SD from the mean), # sex.outlier.sd default = 3
(6) samples with a genotype mismatch (samples with genotype concordance < 0.8 were removed). # sample.genotype.concordance.threshold default = 0.9

# set QC parameters...
qc.parameters <- meffil.qc.parameters(
	beadnum.samples.threshold             = 0.1,
	detectionp.samples.threshold          = 0.1,
	detectionp.cpgs.threshold             = 0.1, 
	beadnum.cpgs.threshold                = 0.1,
	sex.outlier.sd                        = 5,
	snp.concordance.threshold             = 0.95,
	sample.genotype.concordance.threshold = 0.8
)



Probe-wise, we excluded CpGs 
within single nucleotide polymorphisms with a minor allele frequency of >1% located within five nucleotides of the target sites, 
probes tagging non-unique 3΄-subsequences of 30 or more bases long, 
cross-reactive probes as previously described (Chen et al., 2013). 

Additional EPIC-specific cross-reactive probes were obtained from the Maxprobes R package (version 0.0.2) and removed (McCartney et al., 2016; Pidsley et al., 2016) (see Supplementary Figure 4 for an overview of the QC pipeline). 

Functional normalization was performed using the Meffil R package (version 1.3.3). 
Principal component analysis was performed on the control matrix and the 20,000 most variable probes. 
Principal components of the most variable normalized betas corresponding to the iPSC dataset showed a strong batch effect corresponding to the experimental plate.
The ComBat R package was applied only in this dataset to adjust for this batch effect.
