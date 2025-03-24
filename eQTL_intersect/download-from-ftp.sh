#!/usr/bin/env bash
#SBATCH --time 8:00:00
#SBATCH --partition norm

# PhLiPS naive iPSC gene counts expression eQTL
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000023/QTD000399/QTD000399.all.tsv.gz	

# HipSci naive iPSC gene counts expression eQTL
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000016/QTD000361/QTD000361.all.tsv.gz	

#iPSCORE naive iPSC gene counts expression eQTL
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000017/QTD000366/QTD000366.all.tsv.gz	

