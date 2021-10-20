#!/usr/bin/env Rscript

library(UpSetR)

PROJDIR = '/home/tiffinp/epste051/project/select_reseq'

genes = read.table(file.path(PROJDIR, 'results/summary_analyses_and_data_collections/analyses/spanfran2_gwas/gwas_analysis_2021-07-30/candidate_lists/medicago_candidate_windows.tsv'),
                             sep='\t', as.is=TRUE)

phenotypes = c('PC1', 'PC2', 'shannon', 'nod_number',
               'nod_area')

gene_list = unique(genes[genes[, 4] == 'candidate', 5])

data = data.frame('gene'=gene_list,
                  'PC1'=numeric(length(gene_list)),
                  'PC2'=numeric(length(gene_list)),
                  'shannon'=numeric(length(gene_list)),
                  'nod_number'=numeric(length(gene_list)),
                  'nod_area'=numeric(length(gene_list)),
                  stringsAsFactors=FALSE)
for(ph in phenotypes) {
    cands = genes[genes[, 1] == ph & genes[, 4] == 'candidate', 5]
    data[, ph] = as.numeric(data[, 'gene'] %in% cands)
}


overlaps = data[, phenotypes]

pdf('overlap_plot.2021-09-27.pdf', width=12, height=6)
xx = upset(overlaps,
           nsets=ncol(overlaps),
           keep.order=TRUE,
           order.by=c('degree', 'freq'),
           decreasing=c(FALSE, TRUE),
           nintersects=NA,
           empty.intersections=NULL,
           text.scale=1.5,
           set_size.show=TRUE)
print(xx)
dev.off()
