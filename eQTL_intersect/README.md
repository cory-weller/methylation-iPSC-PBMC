
# eQTL Catalog for iPSCs
[https://www.ebi.ac.uk/eqtl/](https://www.ebi.ac.uk/eqtl/)

https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/tabix_ftp_paths.tsv

```bash
# List of ftp filepaths
wget https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv

# PhLiPS naive iPSC gene counts expression eQTL
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000023/QTD000399/QTD000399.all.tsv.gz	

# HipSci naive iPSC gene counts expression eQTL
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000016/QTD000361/QTD000361.all.tsv.gz	

#iPSCORE naive iPSC gene counts expression eQTL
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000017/QTD000366/QTD000366.all.tsv.gz	
```


# Subset to Bonferroni-adjusted P-value significance
```bash
for study in 'QTD000361' 'QTD000366' 'QTD000399'; do
    echo $study
    zcat ${study}.all.tsv.gz | wc -l
    zcat ${study}.all.tsv.gz | awk '($9 < 0.05)' > ${study}.01.tsv  # Temp file to parse later
done
```

|   Study   |  N rows   | Bonferroni-adjusted threshold   |
|-----------|-----------|-----------|
| QTD000361 | 143317269 | `0.05/143317269` = `3.488763e-10` |
| QTD000366 | 142812959 | `0.05/142812959` = `3.501083e-10` |
| QTD000399 | 162075245 | `0.05/162075245` = `3.084987e-10` |

```bash
Rscript filter_threshold.R QTD000361.05.tsv 3.488763e-10
Rscript filter_threshold.R QTD000366.05.tsv 3.501083e-10
Rscript filter_threshold.R QTD000399.05.tsv 3.084987e-10
```

Now the files are small enough to manage in `R`.

