#!/usr/bin/env Rscript
#
# Run PCA with centering but not scaling on median fitness matrices.
# Report both the scores (the "x" component of the prcomp output) and
# the loadings ("rotation" component).
#
# Limit to just the truncatula strains.
#


PROJDIR = paste0(Sys.getenv('HOME'), '/project/select_reseq')
OUTDIR = file.path(PROJDIR, 'results/phenotypes/gwas_phenotypes/spanfran2/2020-11-13_pca')

infile = file.path(PROJDIR, 'results/relative_fitness/compile/2020-05-23/SP2-gwas.tsv')

## Copy script to output directory
cargs = commandArgs(trailingOnly=FALSE)
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUTDIR, 'script_copies'), recursive=TRUE, showWarnings=FALSE)
file.copy(gsub('--file=', '', cargs[grepl('--file=', cargs)]),
          file.path(OUTDIR, 'script_copies'))

truncatula = scan(file.path(PROJDIR, 'data/strain/lists/spanfran2_medicago_truncatula.txt'),
                  what='character')

X = read.csv(infile, sep='\t', header=TRUE, as.is=TRUE, check.names=FALSE)

## Median fitness
med_fit = as.matrix(X[, grepl('_Fit_med', colnames(X))])
colnames(med_fit) = gsub('SP2_', '', gsub('_Fit_med', '', colnames(med_fit)))
rownames(med_fit) = X[, 'strain']
## Added 2021-07-02:
med_fit = med_fit[, truncatula]

## PCA
host_pca = prcomp(t(med_fit), scale=FALSE, center=TRUE)
symbiont_pca = prcomp(med_fit, scale=FALSE, center=TRUE)

## Write output
colnames(symbiont_pca[['x']]) = paste0('rhizobia_score_', colnames(symbiont_pca[['x']]))
write.table(cbind('strain'=rownames(symbiont_pca[['x']]), symbiont_pca[['x']]),
            file=file.path(OUTDIR, 'strain_pca_scores.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
colnames(host_pca[['x']]) = paste0('host_score_', colnames(host_pca[['x']]))
write.table(cbind('strain'=rownames(host_pca[['x']]), host_pca[['x']]),
            file=file.path(OUTDIR, 'plant_pca_scores.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

colnames(symbiont_pca[['rotation']]) = paste0('rhizobia_loading_', colnames(symbiont_pca[['rotation']]))
write.table(cbind('strain'=rownames(symbiont_pca[['rotation']]), symbiont_pca[['rotation']]),
            file=file.path(OUTDIR, 'strain_pca_loadings.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
colnames(host_pca[['rotation']]) = paste0('host_loading_', colnames(host_pca[['rotation']]))
write.table(cbind('strain'=rownames(host_pca[['rotation']]), host_pca[['rotation']]),
            file=file.path(OUTDIR, 'plant_pca_loadings.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

handle = file(file.path(OUTDIR, 'pve.txt'), 'w')
cat('Strain PCA\n',
    paste(symbiont_pca[['sdev']]^2 / sum(symbiont_pca[['sdev']]^2), collapse=', '),
    '\nHost PCA\n',   
    paste(host_pca[['sdev']]^2 / sum(host_pca[['sdev']]^2), collapse=', '),
    sep='', file=handle)
close(handle)

pdf(file.path(OUTDIR, 'host_pca.pdf'), width=7, height=7)
pve = host_pca[['sdev']]^2 / sum(host_pca[['sdev']]^2)
for(i in 1:3) {
    a = (i-1)*2+1
    b = (i-1)*2+2
    plot(host_pca[['x']][, a],
         host_pca[['x']][, b],
         xlab=paste0('PC', a, ' (', round(pve[a],3)*100, '%)'),
         ylab=paste0('PC', b, ' (', round(pve[b],3)*100, '%)'),
         main='Host PCs',
         cex.axis=1.5, cex.lab=1.5)
}
dev.off()

pdf(file.path(OUTDIR, 'rhizobia_pca.pdf'), width=7, height=7)
pve = symbiont_pca[['sdev']]^2 / sum(symbiont_pca[['sdev']]^2)
for(i in 1:3) {
    a = (i-1)*2+1
    b = (i-1)*2+2
    plot(symbiont_pca[['x']][, a],
         symbiont_pca[['x']][, b],
         xlab=paste0('PC', a, ' (', round(pve[a],3)*100, '%)'),
         ylab=paste0('PC', b, ' (', round(pve[b],3)*100, '%)'),
         main='Rhizobia PCs',
         cex.axis=1.5, cex.lab=1.5)
}
dev.off()
