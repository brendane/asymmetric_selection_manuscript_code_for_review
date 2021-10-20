#!/usr/bin/env Rscript
#
# Boxplots of candidate gene values and background values for several
# diversity and divergence statistics.
#

dens_jitter = function(y, x, amount=0.5) {
    h = hist(x, plot=FALSE)
    b = h[['breaks']]
    b[1] = b[1] - 0.01 * diff(range(b))
    splity = split(y, cut(x, breaks=b))
    spliti = split(1:length(y), cut(x, breaks=b))
    result = numeric(length(y))
    #unlist(res_)[unlist(spliti)]
    for(i in seq_along(splity)) {
        xx = splity[[i]]
        d = h[['density']][i] / sum(h[['density']])
        j = jitter(xx, amount=d / 2 * amount)
        result[spliti[[i]]] = j
    }
    result
}

data_dir = '/home/tiffinp/epste051/project/select_reseq/results/summary_analyses_and_data_collections/analyses/spanfran2_gwas/gwas_analysis_2021-07-30'

threshold = 20
stats = c('dvst_beta_max', 'dvst_raisd_max')
stat_names = c('dvst_beta_max'='Max Beta', 'dvst_raisd_max'='Max RAiSD mu')

cols = c('Annot.'='#018571', 'PC1'='#a6611a', 'PC2'='#dfc27d',
         'IG'='#f5f5f5', 'BG'='gray')

data = read.csv(file.path(data_dir, 'rhizobia_results_spreadsheets/rhizobia_data_compiled_by_gene.tsv'),
                sep='\t', header=TRUE, as.is=TRUE, check.names=FALSE, comment.char='#')
genes = rbind(read.csv(file.path(data_dir, 'candidate_lists/gene_list.tsv'),
                       sep='\t', header=FALSE, as.is=TRUE, check.names=FALSE, comment.char='#'),
              read.csv(file.path(data_dir, 'candidate_lists/gene_list_annotated.tsv'),
                       sep='\t', header=FALSE, as.is=TRUE, check.names=FALSE, comment.char='#'))
genes[, 1] = gsub('rhizobia_score_', '', genes[, 1])

perm_data = read.csv(file.path(data_dir, 'comparisons/permutation_comparisons.tsv'),
                     sep='\t', comment.char='#')

cands = function(d, g, p, r, n=FALSE) {
    if(p == 'annotated') {
        result = data[, 'annotated_symbiosis'] == 1 & data[, 'replicon'] == r & data[, 'in_gwas_small'] == 1
    } else {
        if(n) {
            result = data[, 'gene_ID'] %in% g[g[, 1] == p & g[, 3] == r & g[, 4] == 'noncandidate', 5]
        } else {
            result = data[, 'gene_ID'] %in% g[g[, 1] == p & g[, 3] == r & g[, 4] == 'candidate', 5]
        }
    }
    result
}

categories = list(
                  'PC1-chr'=cands(data, genes, 'PC1', 'chromosome'),
                  'PC2-chr'=cands(data, genes, 'PC2', 'chromosome'),
                  'Annot.-chr'=cands(data, genes, 'annotated', 'chromosome'),
                  'BG-chr'=cands(data, genes, 'PC1', 'chromosome', TRUE),
                  'PC1-pSymA'=cands(data, genes, 'PC1', 'psyma'),
                  'PC2-pSymA'=cands(data, genes, 'PC2', 'psyma'),
                  'Annot.-pSymA'=cands(data, genes, 'annotated', 'psyma'),
                  'BG-pSymA'=cands(data, genes, 'PC1', 'psyma', TRUE),
                  'PC1-pSymB'=cands(data, genes, 'PC1', 'psymb'),
                  'PC2-pSymB'=cands(data, genes, 'PC2', 'psymb'),
                  'Annot.-pSymB'=cands(data, genes, 'annotated', 'psymb'),
                  'BG-pSymB'=cands(data, genes, 'PC1', 'psymb', TRUE)
                  )


pdf('boxplot_popgen_stats.all_replicons.2021-10-14.pdf', width=7, height=8.5)
par(mfcol=c(3, 2), mar=par('mar')/1.2, mgp=par('mgp')/1.3, bty='l')

for(stat in stats) {

    for(replicon in c('chr', 'pSymA', 'pSymB')) {

        repl_categories = categories[grep(replicon, names(categories))]
        stat_data = structure(vector('list', length=length(repl_categories)),
                              names=names(repl_categories))
        for(ct in names(stat_data)) {
            if(is.character(data[, stat])) {
                s = structure(as.numeric(unlist(strsplit(data[repl_categories[[ct]], stat], ','))),
                              names=data[repl_categories[[ct]], 'gene_ID'])
            } else {
                s = structure(data[repl_categories[[ct]], stat],
                              names=data[repl_categories[[ct]], 'gene_ID'])
            }
            s = s[!is.na(s) & !is.infinite(s)]
            stat_data[[ct]] = s
        }

        bx = boxplot(stat_data,
                     col=ifelse(grepl('BG', names(stat_data)), cols['BG'], NA),
                     border=ifelse(grepl('BG', names(stat_data)), 'black', NA),
                     varwidth=FALSE,
                     outline=TRUE,
                     las=2,
                     cex.axis=1.5,
                     cex.lab=1.5,
                     ylab=stat_names[stat],
                     ylim={ if(replicon == 'chr' && stat == 'dvst_beta_max') c(-2, 11) else NULL },
                     names = gsub('-.+', '', names(stat_data)),
                     xgap.axis=0,
                     ygap.axis=0,
                     xpd=NA,
                     main=replicon
                     )

        for(i in seq_along(stat_data)) {
            ph = gsub('-.+', '', names(stat_data)[i])
            if(grepl('BG', ph)) next
            cl = cols[ph]
            if(ph == 'IG') cl = 'black'
            x = stat_data[[i]]
            points(dens_jitter(rep(i, length(x)), x, amount=0.5), x,
                   col=cl, cex=1.2, pch=1, xpd=NA)
            segments(i-0.5, median(x), i+0.5, median(x), lwd=4, col=cl)
            cat('median of', stat, replicon, ph, median(x), '\n')
        }

        ## Mean of PC1 and PC2 permutation medians (they are very similar
        ## to each other)
        md = mean(perm_data[perm_data[, 'statistic'] == stat &
                            perm_data[, 'phenotype'] %in% c('PC1', 'PC2') &
                            perm_data[, 'replicon'] == tolower(gsub('chr', 'chromosome', replicon)),
                            'perm_median'])
                       
        abline(h=md, col='gray30', lty=2, lwd=3)
        print(c(stat, replicon))
        print(sort(stat_data[[4]], decreasing=TRUE)[1:10])
    }
}

dev.off()

