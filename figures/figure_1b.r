#!/usr/bin/env Rscript
#
# Distribution of median relative fitness
#
# UPDATE 23 August 2021: heatmaps ordered by PC scores and loadings
#

library(heatmap3)
library(magrittr)
library(psych)

to_named_vector = function(x, a, b) {
    structure(x[, a], names=x[, b])
}

rpcs = read.csv('/home/tiffinp/epste051/project/select_reseq/results/phenotypes/gwas_phenotypes/spanfran2/2020-11-13_pca/strain_pca_scores.tsv',
               sep='\t')
rpc_loadings = read.csv('/home/tiffinp/epste051/project/select_reseq/results/phenotypes/gwas_phenotypes/spanfran2/2020-11-13_pca/strain_pca_loadings.tsv',
                        sep='\t')
rpc1 = structure(rpcs[, 'rhizobia_score_PC1'], names=rpcs[, 'strain'])
rpc1l = structure(rpc_loadings[, 'rhizobia_loading_PC1'], names=rpc_loadings[, 'strain'])
rpc2 = structure(rpcs[, 'rhizobia_score_PC2'], names=rpcs[, 'strain'])
rpc2l = structure(rpc_loadings[, 'rhizobia_loading_PC2'], names=rpc_loadings[, 'strain'])
mpcs = read.csv('/home/tiffinp/epste051/project/select_reseq/results/phenotypes/gwas_phenotypes/spanfran2/2020-11-13_pca/plant_pca_scores.tsv',
               sep='\t')
mpc_loadings = read.csv('/home/tiffinp/epste051/project/select_reseq/results/phenotypes/gwas_phenotypes/spanfran2/2020-11-13_pca/plant_pca_loadings.tsv',
                        sep='\t')
mpc1 = structure(mpcs[, 'host_score_PC1'], names=mpcs[, 'strain'])
mpc1l = structure(mpc_loadings[, 'host_loading_PC1'], names=mpc_loadings[, 'strain'])
mpc2 = structure(mpcs[, 'host_score_PC2'], names=mpcs[, 'strain'])
mpc2l = structure(mpc_loadings[, 'host_loading_PC2'], names=mpc_loadings[, 'strain'])

truncatula = scan('/home/tiffinp/epste051/project/select_reseq/data/strain/lists/spanfran2_medicago_truncatula.txt',
                  what='character')
shannon = read.csv('~/project/select_reseq/results/phenotypes/gwas_phenotypes/spanfran2/2021-06-29_nod_shannon/phenotypes.tsv',
                   sep='\t', as.is=TRUE) %>%
    to_named_vector('shannon', 'strain') %>%
    .[truncatula] %>%
    sort()

fit_data = read.csv('/home/tiffinp/epste051/project/select_reseq/results/relative_fitness/spanfran2/2020-10-27/real/phenotypes.tsv',
                    sep='\t', header=TRUE, as.is=TRUE, check.names=FALSE)
rel_fit = as.matrix(fit_data[, grepl('_Fit_med', colnames(fit_data))])
colnames(rel_fit) = gsub('_Fit_med', '', colnames(rel_fit))
rownames(rel_fit) = fit_data[, 'strain']
rel_fit = rel_fit[, truncatula]

mean_rel_fit = rowMeans(rel_fit)
max_rel_fit = apply(rel_fit, 1, max)
var_rel_fit = apply(rel_fit, 2, var)

## Because of negative numbers, this is calculated
## through some transformations
hmean_rel_fit = apply(rel_fit, 1, function(x) log2(harmonic.mean(2^x)))

pdf('rhizobial_relative_fitness_distribution.2021-07-22.pdf', width=3.7, height=9)
par(mar=c(2, 4, 2, 2)+0.1, mgp=par('mgp')/1.5)
layout(rbind(1, 2, 2, 2))

mx = 4.5
mn = -8

hist(mean_rel_fit, breaks=seq(mn, mx, 0.5), col='gray20', border='white',
     xaxs='i', yaxs='i', cex.axis=1.5, cex.lab=1.5, main='',
     ylab='Number of strains',
     xlab='Mean of median relative fitness\nacross host genotypes',
     las=1)

plot(0, 0, xlim=c(mn, mx), ylim=c(0, nrow(rel_fit)+1),
    type='n', xlab='Median relative fitness', xaxs='i', xaxt='n',
    ylab='', yaxs='i', yaxt='n', bty='n',
    cex.axis=1.5, cex.lab=1.5)
axis(side=1, lwd=1, lwd.ticks=1, cex.axis=1.5)

sorted_rel_fit = rel_fit[order(mean_rel_fit), ]
for(i in 1:nrow(sorted_rel_fit)) {
    #if(i %% 2 == 0) {
    #    rect(mn+0.05, i-0.4, mx, i+0.4, col=rgb(0.99, 0.96, 0.94), border=NA)
    #}
    colr = c('#b360d14D', '#6eb8304D')[i %% 2 + 1]
           #col=rgb(0.4, 0.6, 0.9, 0.1),
    m = min(sorted_rel_fit[i, ])
    #segments(mn, i, m, i, col='gray70', lwd=0.5)
    points(sorted_rel_fit[i, ], jitter(rep(i, ncol(sorted_rel_fit)), amount=0.2),
           col=colr,
           lwd=0,
           cex=0.5, pch=19, xpd=NA)
    #text({if(m > -5) m else max(sorted_rel_fit[i, ])},
    #    i, rownames(sorted_rel_fit)[i], cex=0.65,
    #     col=c('#8100e6', '#2b9900')[i %% 2 + 1],
    #     pos={if(m > -5) 2 else 4}, offset=1)
    axis(side=2, lwd=0, lwd.ticks=0, gap.axis=0,
         cex.axis=0.65, at=i, col.ticks='gray70',
         labels=rownames(sorted_rel_fit)[i], col.axis=c('#8100e6', '#2b9900')[i %% 2 + 1], las=1)
} 
abline(v=0, lty=2, col='gray20')
#axis(side=2, lwd=0, lwd.ticks=1, gap.axis=0,
#     cex.axis=0.7, at=1:nrow(sorted_rel_fit),
#     labels=rownames(sorted_rel_fit), las=1)

dev.off()




pdf('medicago_relative_fitness_distribution.2021-07-22.pdf', width=3.7, height=9)

double_sorted_rel_fit = t(rel_fit[order(mean_rel_fit),
                                  match(names(shannon), colnames(rel_fit))])

m = colMeans(double_sorted_rel_fit)
m = (m - min(m)) / (max(m) - min(m))
pal = c("#3a4d6b", "#3d6da2", "#799a96", "#ccbe6a", "#ffec99")
cr = colorRamp(pal)
#cr = colorRamp(c('blue', 'blue'))
cl = apply(cr(m), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255, 0.2))

pm = par('mar')
par(mar=c(pm[1], pm[2], 0.5, pm[4]))

plot(0, 0, xlim=c(mn, mx), ylim=c(0, nrow(double_sorted_rel_fit)+1),
    type='n', xlab='Median relative fitness', xaxs='i', xaxt='n',
    ylab='', yaxs='i', yaxt='n', bty='l',
    cex.axis=1.5, cex.lab=1.5)
axis(side=1, lwd=0, lwd.ticks=1, cex.axis=1.5, gap.axis=0)

for(i in 1:nrow(double_sorted_rel_fit)) {
    if(i %% 2 == 0) {
        rect(mn+0.05, i-0.4, mx, i+0.4, col='gray95', border=NA)
    }
    
    points(double_sorted_rel_fit[i, ], rep(i, ncol(double_sorted_rel_fit)),
           col=cl,
           lwd=0,
           cex=0.35, pch=19)
} 
axis(side=2, lwd=0, lwd.ticks=1, gap.axis=0, cex.lab=0.25,
     cex.axis=0.25, at=1:nrow(double_sorted_rel_fit),
     labels=rownames(double_sorted_rel_fit), las=1)

dev.off()


pdf('rhizobial_fitness_tradeoffs_2021-07-22.pdf', width=6, height=3)

par(mfrow=c(1, 2))
n_gt_0 = apply(rel_fit, 1, function(x) sum(x > 0))
plot(n_gt_0, max_rel_fit, col='gray20',
     cex.axis=1.5, cex.lab=1.5, main='', xlab='# Hosts fitness > 0',
     ylab='Max relative fitness', las=1, cex=0.7, xpd=NA, bty='L')
legend('topleft', fill=NA, legend=paste('r =', round(cor(n_gt_0, max_rel_fit), 2)),
       border=NA, bty='n')
plot(hmean_rel_fit, max_rel_fit, col='gray20',
     cex.axis=1.5, cex.lab=1.5, main='', xlab='Harmonic mean relative\nfitness',
     ylab='', las=1, cex=0.7, xpd=NA, bty='L')
legend('topleft', fill=NA, legend=paste('r =', round(cor(hmean_rel_fit, max_rel_fit), 2)),
       border=NA, bty='n')

dev.off()

pdf('fitness_heatmap_2021-07-22.pdf', width=6, height=6)

stopifnot(all(names(rpc1) == rownames(rel_fit)))
stopifnot(all(names(rpc2) == rownames(rel_fit)))
stopifnot(all(names(mpc1) == colnames(rel_fit)))
stopifnot(all(names(mpc2) == colnames(rel_fit)))

rcc1 = apply(colorRamp(c('red', 'gray', 'blue'))(rank(rpc1)/length(rpc1)),
           1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
rcc2 = apply(colorRamp(c('red', 'gray', 'blue'))(rank(rpc2)/length(rpc2)),
           1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
mcc1 = apply(colorRamp(c('red', 'gray', 'blue'))(rank(mpc1)/length(mpc1)),
           1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
mcc2 = apply(colorRamp(c('red', 'gray', 'blue'))(rank(mpc2)/length(mpc2)),
           1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))


heatmap3(rel_fit, scale='none', distfun=function(x) dist(x, method='euclidean'),
        col=colorRampPalette(pal)(100),
        cexRow=0.25, cexCol=0.25, balanceColor=TRUE)

heatmap3(t(apply(rel_fit, 1, function(x) x-mean(x))),
         scale='none', distfun=function(x) dist(x, method='euclidean'),
         col=colorRampPalette(pal)(100),
         cexRow=0.25, cexCol=0.25, balanceColor=TRUE)

dev.off()

centerr = function(x) t(apply(x, 1, function(y) y - mean(y)))
centerc = function(x) apply(x, 2, function(y) y - mean(y))

pdf('fitness_heatmap_ordered_by_PCs.pdf', width=10, height=10)

heatmap3(rel_fit[names(sort(mpc1l)), names(sort(mpc1))],
         scale='none', Rowv=NA, Colv=NA,
         col=colorRampPalette(pal)(100),
         cexRow=0.25, cexCol=0.25, balanceColor=TRUE,
         main='Host PC1')

heatmap3(rel_fit[names(sort(mpc2l)), names(sort(mpc2))],
         scale='none', Rowv=NA, Colv=NA,
         col=colorRampPalette(pal)(100),
         cexRow=0.25, cexCol=0.25, balanceColor=TRUE,
         main='Host PC2')

heatmap3(rel_fit[names(sort(rpc1)), names(sort(rpc1l))],
         scale='none', Rowv=NA, Colv=NA,
         col=colorRampPalette(pal)(100),
         cexRow=0.25, cexCol=0.25, balanceColor=TRUE,
         main='Rhizobia PC1')

heatmap3(rel_fit[names(sort(rpc2)), names(sort(rpc2l))],
         scale='none', Rowv=NA, Colv=NA,
         col=colorRampPalette(pal)(100),
         cexRow=0.25, cexCol=0.25, balanceColor=TRUE,
         main='Rhizobia PC2')

dev.off()

#heatmap3(rel_fit, scale='none', distfun=function(x) dist(x, method='euclidean'),
#        col=colorRampPalette(pal)(100),
#        RowSideColors=rcc1, ColSideColors=mcc1,
#        cexRow=0.25, cexCol=0.25, balanceColor=TRUE)
#
#heatmap3(rel_fit, scale='none', distfun=function(x) dist(x, method='euclidean'),
#        col=colorRampPalette(pal)(100),
#        RowSideColors=rcc2, ColSideColors=mcc2,
#        cexRow=0.25, cexCol=0.25, balanceColor=TRUE)
