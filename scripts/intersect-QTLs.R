#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)

eqtl_files <- c("eQTL_intersect/QTD000361.all.tsv.gz",
                "eQTL_intersect/QTD000366.all.tsv.gz",
                "eQTL_intersect/QTD000399.all.tsv.gz"
)

eqtl <- fread(eqtl_files[1])