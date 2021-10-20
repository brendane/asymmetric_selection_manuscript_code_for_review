#!/usr/bin/env Rscript
#
# Boxplots of candidate window values and background values for several
# diversity and divergence statistics.
#

data_dir = '/home/tiffinp/epste051/project/select_reseq/results/summary_analyses_and_data_collections/analyses/spanfran2_gwas/gwas_analysis_2021-07-30'

threshold = 0.01
stats = c('max_beta_selection', 'max_raisd_mu_score')
stat_names = c('max_beta_selection'='Max Beta', 'max_raisd_mu_score'='Max RAiSD mu')

cols = c('Annot.'='#a40122',
         'PC1'='#9f0162',
         'PC2'='#00c2f9',
         'Shannon'='#e20134',
         'Nod Number'='#00fccf',
         'Nod Area'='#8400cd',
         'BG'='gray')

#'Nod Mass'='#ff5aaf',

rank_columns = c('PC1'='host_score_PC1_max_effect_size_lmm',
                 'PC2'='host_score_PC2_max_effect_size_lmm',
                 'Shannon'='shannon_max_effect_size_lmm',
                 'Nod Number'='nod_number_max_effect_size_lmm',
                 'Nod Area'='nod_area_max_effect_size_lmm',
                 'Annot.'='annotated_symbiosis')

data = read.csv(file.path(data_dir, 'medicago_results_spreadsheets/medicago_data_compiled_by_window.tsv'),
                sep='\t', header=TRUE, as.is=TRUE, check.names=FALSE, comment.char='#')
rownames(data) = paste(data[, 'chromosome'], data[, 'start'], data[, 'end'], sep='-')


candidates = structure(vector('list', length=length(rank_columns)+1),
                       names=c(names(rank_columns), 'BG'))
cand_data = read.csv(file.path(data_dir, 'candidate_lists/medicago_candidate_windows.tsv'),
                     header=FALSE, as.is=TRUE, sep='\t')
for(pheno in names(candidates)) {
    if(pheno == 'BG') {
        candidates[[pheno]] = unique(cand_data[, 5])
    } else if(pheno == 'Annot.') {
        candidates[[pheno]] = unique(cand_data[cand_data[, 1] == 'annotated' &
                                               cand_data[, 4] == 'candidate',
                                               5])
    } else {
        if(grepl('^PC', pheno)) {
            p = pheno
        } else {
            p = gsub(' ', '_', tolower(pheno))
        }
        candidates[[pheno]] = cand_data[cand_data[, 1] == p &
                                        cand_data[, 2] == 'effect_size' &
                                        cand_data[, 3] == 0.01 &
                                        cand_data[, 4] == 'candidate',
                                        5]
    }
}


pdf('medicago_boxplot_popgen_stats.2021-09-27.pdf', width=9, height=6)

par(mfcol=c(1, 2), mgp=par('mgp')/1.3, bty='l')

for(stat in stats) {

    stat_data = structure(vector('list', length=length(candidates)),
                          names=names(candidates))
    max_wh = 0
    for(ct in names(stat_data)) {
        s = data[candidates[[ct]], stat]
        s = s[!is.na(s) & !is.infinite(s)]
        s[s < 0] = 0
        iqr = quantile(s, c(0.25, 0.75))
        wh = diff(iqr) * 1.5 + iqr[2]
        stat_data[[ct]] = s
        if(wh > max_wh) max_wh = wh
    }


    for(ct in names(stat_data)) {
        d = stat_data[[ct]]
        dd = d[d > max_wh]
        if(length(dd) > 0) {
            d[d > max_wh] = max_wh + log2(dd - max_wh)
            stat_data[[ct]] = d
        }
    }

    bx = boxplot(stat_data,
                 col=cols[names(stat_data)],
                 outcol=cols[names(stat_data)],
                 border='black',
                 varwidth=FALSE,
                 outline=TRUE,
                 las=2,
                 cex.axis=1.5,
                 cex.lab=1.5,
                 cex=0.5,
                 ylab=stat_names[stat],
                 names = gsub('-.+', '', names(stat_data)),
                 xgap.axis=0,
                 ygap.axis=0,
                 xpd=NA,
                 yaxt='n'
    )

    stp = 20
    if(max(unlist(data[, stat]), na.rm=TRUE) > 100) {
        stp = 50
    }
    at = c(seq(-5, max_wh, 5), seq(stp, max(unlist(data[, stat]), na.rm=TRUE), stp))
    atl = at
    at[at > max_wh] = max_wh + log2(at[at > max_wh] - max_wh)

    axis(side=2, at=at, labels=atl)
    points(par('usr')[1], max_wh, xpd=NA)

    print(lapply(stat_data, function(x) sum(!is.na(x) & !is.infinite(x))))

    abline(h=median(stat_data[['BG']]), lwd=3, lty=2, col='gray30')

}

dev.off()
