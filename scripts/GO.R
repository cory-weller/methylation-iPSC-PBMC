#!/usr/bin/env Rscript


package_list = c('ggplot2', 'data.table', 'corrplot', 'umap', 'magick', 'ggdendro', 'ecodist','ggbeeswarm', 'ggrepel', 'ggthemes', 'foreach','reshape2','org.Hs.eg.db','clusterProfiler','pheatmap')
lapply(package_list, require, character.only=TRUE)



enrich_pvalue <- 0.05

all_gene_vector <- unique(as.character(dat$lead_gene))

dat <- fread('methQTL/PBMC.cis_qtl.txt.gz')

enriched_set <- 
if (!dir.exists(outdir)){
    dir.create(outdir,recursive = T)
}

    #if (nrow(kk)>0){
    #  barplot(kk, showCategory=20)
    #  ggsave(file.path(outdir,paste0(out_prefix,'.kegg.pdf')),width=width,height=height)
    # }
    print('Save enrichment analysis results')
    ezwrite(enrich_res, outdir, paste0(out_prefix,'.enrich_res.tsv'))
    #write.csv(enrich_res,file.path(outdir,paste0(out_prefix,'.enrich_res.csv')),row.names = F)