#!/usr/bin/env Rscript
#
# Make Venn Diagrams showing overlap among Medicago candidates:
#   - different phenotypes
#   - LM and LMM analyses
#
# The units of the analysis are windows.
#

phenotypes = c('PC1', 'PC2')
replicons = c('chromosome', 'psyma', 'psymb')

color = c('PC1'='red', 'PC2'='blue')

data_dir = '/home/tiffinp/epste051/project/select_reseq/results/summary_analyses_and_data_collections/analyses/spanfran2_gwas/gwas_analysis_2021-07-30'

cand_genes = read.csv(file.path(data_dir, 'candidate_lists/list_of_rhizobia_candidates.tsv'),
                     sep='\t', header=FALSE, as.is=TRUE)
gene_data = read.csv(file.path(data_dir, 'rhizobia_results_spreadsheets/rhizobia_data_compiled_by_gene.tsv'),
                     comment.char='#', sep='\t')

combns = apply(expand.grid(phenotypes, replicons), 1, paste, collapse='-')

candidates = structure(vector('list', length=length(combns)), names=combns)
lm_candidates = candidates
for(cb in names(candidates)) {
    pheno = unlist(strsplit(cb, '-'))[1]
    replicon = unlist(strsplit(cb, '-'))[2]
    candidates[[cb]] = cand_genes[cand_genes[, 1] == paste0('rhizobia_score_', pheno) &
                                  cand_genes[, 2] == 'effect_size' &
                                  cand_genes[, 3] == replicon &
                                  cand_genes[, 4] == 'candidate', 5]
    lm_candidates[[cb]] = cand_genes[cand_genes[, 1] == paste0('rhizobia_score_', pheno, '_lm') &
                                     cand_genes[, 2] == 'effect_size' &
                                     cand_genes[, 3] == replicon &
                                     cand_genes[, 4] == 'candidate', 5]
}


do_venn_diagram_plot = function(candidates1, candidates2) {
    phenos = intersect(names(candidates1), names(candidates2))
    np = length(phenos)

    bars = structure(numeric(length(phenos)), names=phenos)
    for(ph in phenos) {
        bars[ph] = length(intersect(candidates1[[ph]], candidates2[[ph]]))
    }

    barplot(bars, col=color[sapply(names(bars), function(x) unlist(strsplit(x, '-'))[1])],
            border='white',
            xlab='', ylab='Number of candidate windows\nin both LM and LMM',
            gap.axis=0, las=1)
}


### Plot number of candidates in both LM and LMM
#pdf('ensifer_candidates_lm_lmm_venn_diagram.2021-10-01.pdf', width=6, height=6)
#do_venn_diagram_plot(candidates, lm_candidates)
#dev.off()


## Plot number of candidates in LMM as minimum LM ranking
## increases (up to 50 genes)

pdf('ensifer_lmm_candidates_in_more_than_20_lm.pdf', width=6, height=10)
par(mfrow=c(3, 2))
for(cb in names(candidates)) {
    pheno = unlist(strsplit(cb, '-'))[1]
    replicon = unlist(strsplit(cb, '-'))[2]
    cn = paste0('max_effect_size_small_groups_', pheno, '_lm_', gsub('omosome', '', replicon))
    lm_ranked = gene_data[order(gene_data[, cn], decreasing=TRUE, na.last=TRUE), 'gene_ID']
    nc = length(candidates[[cb]])
    n_lmm_cands = numeric(50)
    n = 0
    for(i in 1:length(n_lmm_cands)) {
        if(lm_ranked[i] %in% candidates[[cb]]) n = n+1
        n_lmm_cands[i] = n
    }

    plot(n_lmm_cands, xlab='LM candidate cutoff',
         ylab='Number of LMM candidates included',
         col=color[pheno],
         ylim=c(0, nc),
         main=gsub('-', ': ', cb),
         cex.axis=1.5, cex.lab=1.5, las=1,
         xaxs='i', yaxs='i')
}
dev.off()
