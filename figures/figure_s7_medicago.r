#!/usr/bin/env Rscript
#
# Make Venn Diagrams showing overlap among Medicago candidates:
#   - different phenotypes
#   - LM and LMM analyses
#
# The units of the analysis are windows.
#

phenotypes = c('PC1', 'PC2', 'shannon', 'nod_number', 'nod_area')
threshold = 0.01

color = c('PC1'='#19aeff', 'PC2'='#005c94', 'shannon'='#b88100',
          'nod_number'='#b50000', 'nod_area'='#ff6600')

data_dir = '/home/tiffinp/epste051/project/select_reseq/results/summary_analyses_and_data_collections/analyses/spanfran2_gwas/gwas_analysis_2021-07-30'

cand_wins = read.csv(file.path(data_dir, 'candidate_lists/medicago_candidate_windows.tsv'),
                     sep='\t', header=FALSE, as.is=TRUE)

window_data = read.csv(file.path(data_dir, 'medicago_results_spreadsheets/medicago_data_compiled_by_window.tsv'),
                       sep='\t', comment.char='#')


candidates = structure(vector('list', length=length(phenotypes)),
                       names=phenotypes)
lm_candidates = structure(vector('list', length=length(phenotypes)),
                          names=phenotypes)
for(pheno in names(candidates)) {
    candidates[[pheno]] = cand_wins[cand_wins[, 1] == pheno &
                                    cand_wins[, 2] == 'effect_size' &
                                    cand_wins[, 3] == threshold &
                                    cand_wins[, 4] == 'candidate', 5]
    if(! pheno %in% c('annotated', 'IG')) {
        lm_candidates[[pheno]] = cand_wins[cand_wins[, 1] == paste0(pheno, '_lm') &
                                           cand_wins[, 2] == 'effect_size' &
                                           cand_wins[, 3] == threshold &
                                           cand_wins[, 4] == 'candidate', 5]
    }
}


do_venn_diagram_plot = function(candidates1, candidates2) {
    phenos = intersect(names(candidates1), names(candidates2))
    np = length(phenos)

    bars = structure(numeric(length(phenos)), names=phenos)
    for(ph in phenos) {
        bars[ph] = length(intersect(candidates1[[ph]], candidates2[[ph]]))
    }

    barplot(bars, col=color[names(bars)], border='white',
            xlab='', ylab='Number of candidate windows\nin both LM and LMM',
            gap.axis=0, las=1)
}


## Plot number of candidates in both LM and LMM
pdf('medicago_candidates_lm_lmm_venn_diagram.2021-10-01.pdf', width=6, height=6)
do_venn_diagram_plot(candidates, lm_candidates)
dev.off()

## Plot number of candidates in LMM as minimum LM ranking
## increases (up to 4%)

pdf('medicago_lmm_candidates_in_more_than_1percent_lm.pdf', width=6, height=10)
par(mfrow=c(3, 2))
for(phenotype in phenotypes) {
    cn = paste0(gsub('PC', 'host_score_PC', phenotype), '_max_effect_size_lm')
    lm_ranked = paste(window_data[, 'chromosome'], window_data[, 'start'],
                      window_data[, 'end'], sep='-')[order(window_data[, cn], decreasing=TRUE, na.last=TRUE)]
    nc = length(candidates[[phenotype]])
    n_lmm_cands = numeric(floor(nc*4))
    n = 0
    for(i in 1:length(n_lmm_cands)) {
        if(lm_ranked[i] %in% candidates[[phenotype]]) n = n+1
        n_lmm_cands[i] = n
    }

    plot(n_lmm_cands, xlab='LM candidate cutoff',
         ylab='Number of LMM candidates included',
         col=color[phenotype],
         ylim=c(0, nc),
         main=phenotype,
         cex.axis=1.5, cex.lab=1.5, las=1,
         xaxs='i', yaxs='i',
         xaxt='n')
    axis(side=1, at=seq(nc, nc*4, nc),
         label=c('1%', '2%', '3%', '4%'))
    abline(v=seq(nc, nc*4, nc), lty=3)
    abline(h=nc, lty=3)
}
dev.off()
