#!/usr/bin/env Rscript
#
# PC1 and PC2 scores for plants and rhizobia
#

## TODO: show correlation between average fitness and first PC

mtr = read.csv('/home/tiffinp/epste051/project/select_reseq/results/phenotypes/gwas_phenotypes/spanfran2/2020-11-13_pca/plant_pca_scores.tsv',
               sep='\t', header=TRUE, check.names=FALSE, as.is=TRUE) 
rhiz = read.csv('/home/tiffinp/epste051/project/select_reseq/results/phenotypes/gwas_phenotypes/spanfran2/2020-11-13_pca/strain_pca_scores.tsv',
                sep='\t', header=TRUE, check.names=FALSE, as.is=TRUE) 
pves = scan('/home/tiffinp/epste051/project/select_reseq/results/phenotypes/gwas_phenotypes/spanfran2/2020-11-13_pca/pve.txt',
            what='character', sep='\n')

truncatula = scan('/home/tiffinp/epste051/project/select_reseq/data/strain/lists/spanfran2_medicago_truncatula.txt',
                  what='character')
fit_data = read.csv('/home/tiffinp/epste051/project/select_reseq/results/relative_fitness/spanfran2/2020-10-27/real/phenotypes.tsv',
                    sep='\t', header=TRUE, as.is=TRUE, check.names=FALSE)
rel_fit = as.matrix(fit_data[, grepl('_Fit_med', colnames(fit_data))])
colnames(rel_fit) = gsub('_Fit_med', '', colnames(rel_fit))
rownames(rel_fit) = fit_data[, 'strain']
rel_fit = rel_fit[, truncatula]
median_rel_fit = apply(rel_fit, 1, median)

pve_rhiz = numeric(ncol(rhiz)-1)
pve_mtr = numeric(ncol(mtr)-1)
r = FALSE
m = FALSE
for(l in pves) {
    if(r) {
        pve_rhiz = as.numeric(unlist(strsplit(l, ', ')))
        r = FALSE
    }
    if(m) {
        pve_mtr = as.numeric(unlist(strsplit(l, ', ')))
        m = FALSE
    }
    if(l == 'Strain PCA') {
        r = TRUE
        m = FALSE
    }
    if(l == 'Host PCA') {
        m = TRUE
        r = FALSE
    }
}


pdf('pc_score_plots.2021-07-22.pdf', width=3, height=9)
par(mfrow=c(3, 1), mgp=par('mgp')/1.5)
plot(mtr[, 'host_score_PC1'], mtr[, 'host_score_PC2'],
     cex=1, col='gray20', cex.axis=1.5,
     cex.lab=1.5,
     xlab=paste0('PC1 (', round(pve_mtr[1], 2)*100, '%)'),
     ylab=paste0('PC2 (', round(pve_mtr[2], 2)*100, '%)'),
     main='Host', ygap.axis=0, xgap.axis=0)
plot(rhiz[, 'rhizobia_score_PC1'], rhiz[, 'rhizobia_score_PC2'],
     cex=1, col='gray20', cex.axis=1.5,
     cex.lab=1.5,
     xlab=paste0('PC1 (', round(pve_rhiz[1], 2)*100, '%)'),
     ylab=paste0('PC2 (', round(pve_rhiz[2], 2)*100, '%)'),
     main='Rhizobia', ygap.axis=0, xgap.axis=0)
plot(rhiz[, 'rhizobia_score_PC1'], median_rel_fit[rhiz[, 'strain']],
     cex=1, main='PC1 and Fitness',
     xlab='PC1', ylab='Median fitness across hosts',
     cex.lab=1.5, cex.axis=1.5)
dev.off()


pdf('pc_score_plots.2021-07-22.labeled.pdf', width=3, height=6)
par(mfrow=c(2, 1), mgp=par('mgp')/1.5)
plot(mtr[, 'host_score_PC1'], mtr[, 'host_score_PC2'],
     cex=1, col='gray20', cex.axis=1.5,
     cex.lab=1.5, type='n', bty='L',
     xlab=paste0('PC1 (', round(pve_mtr[1], 2)*100, '%)'),
     ylab=paste0('PC2 (', round(pve_mtr[2], 2)*100, '%)'),
     main='Host', ygap.axis=0, xgap.axis=0)
text(mtr[, 'host_score_PC1'], mtr[, 'host_score_PC2'],
     mtr[, 'strain'], cex=0.3, xpd=NA)
plot(rhiz[, 'rhizobia_score_PC1'], rhiz[, 'rhizobia_score_PC2'],
     cex=1, col='gray20', cex.axis=1.5,
     cex.lab=1.5, type='n', bty='L',
     xlab=paste0('PC1 (', round(pve_rhiz[1], 2)*100, '%)'),
     ylab=paste0('PC2 (', round(pve_rhiz[2], 2)*100, '%)'),
     main='Rhizobia', ygap.axis=0, xgap.axis=0)
text(rhiz[, 'rhizobia_score_PC1'], rhiz[, 'rhizobia_score_PC2'],
     rhiz[, 'strain'], cex=0.3, xpd=NA)
dev.off()
