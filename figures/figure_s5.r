#!/usr/bin/env Rscript
#
# Distribution of LD group size per replicon
#

library(magrittr)

PROJDIR = '/home/tiffinp/epste051/project/'

replicon_colors = c('chromosome'=rgb(0.105882353, 0.6196078, 0.466666667),
                    'psymb'=rgb(0.458823529, 0.4392157, 0.701960784),
                    'psyma'=rgb(0.850980392, 0.3725490, 0.007843137))

ld_groups_file = file.path(PROJDIR, 'select_reseq/results/ld/pangenome/r2_groups/2021-04-06_denovo/SP2/0.95/output.tsv')
ld_group_replicon_file = file.path(PROJDIR, 'select_reseq/results/summary_analyses_and_data_collections/analyses/spanfran2_gwas/gwas_analysis_2021-07-30/ld_group_replicons.tsv')

ld_groups = read.csv(ld_groups_file, sep='\t', as.is=TRUE, header=TRUE)
ld_groups[, 'group'] = paste0('group-', ld_groups[, 'group'])
ld_group_replicons = read.csv(ld_group_replicon_file, sep='\t', header=TRUE, as.is=TRUE)
names(ld_group_replicons)[1] = 'ldg'

ld_group_nvars = aggregate(ld_groups[, 'rs'], list(ld_groups[, 'group']), length)
names(ld_group_nvars) = c('ldg', 'nvars')
ld_group_ngenes = aggregate(ld_groups[, 'rs'], list(ld_groups[, 'group']),
                            function(rs) length(unique( gsub('-.+', '', gsub('snp-|sibeliaz-', '', rs)) )))
names(ld_group_ngenes) = c('ldg', 'ngenes')
ld_group_pav_prop = aggregate(ld_groups[, 'rs'], list(ld_groups[, 'group']),
                              function(rs) sum(grepl('sibeliaz', rs)) / length(rs))
names(ld_group_pav_prop) = c('ldg', 'pavp')

ld_data = merge(ld_group_ngenes, ld_group_pav_prop) %>%
    merge(ld_group_nvars) %>%
    merge(ld_group_replicons)

pdf('group_size_distribution_figure.2021-08-02.pdf', width=6, height=5)
    
par(mfrow=c(3,3), mgp=par('mgp')/1.5, 'mar'=c(2.1, 4.1, 1.5, 1.5))
for(replicon in c('chromosome', 'psyma', 'psymb')) {

    color = replicon_colors[[replicon]]

    h = hist(ld_data[ld_data[, 'assignment'] == replicon, 'nvars'],
             breaks=seq(0, 2000, 1), plot=FALSE)
    h[['counts']] = c(h[['counts']][1:20], sum(tail(h[['counts']], -20)))
    h[['mids']] = h[['mids']][1:21]
    h[['breaks']] = h[['breaks']][1:22]
    plot(h, xlim=c(0, 21),
         col=color, border='white', cex.axis=1.5, cex.lab=1.5,
         xlab='Number of variants', ylab='', main='',
         yaxs='i', xaxs='i', lwd.axis=0, lwd.ticks=1, bty='l',
         xaxt='n', las=1, xpd=NA)
    axis(side=1, at=c(0, 5, 10, 15, 20.5),
         labels=c('0', '5', '10', '15', '>20'),
         cex.axis=1.5, lwd=0, lwd.ticks=1, gap.axis=0, xpd=NA)
    grid(nx=NA, ny=NULL)
    box(bty='l')

    h = hist(ld_data[ld_data[, 'assignment'] == replicon, 'ngenes'],
             breaks=seq(0, 2000, 1), plot=FALSE)
    h[['counts']] = c(h[['counts']][1:20], sum(tail(h[['counts']], -20)))
    h[['mids']] = h[['mids']][1:21]
    h[['breaks']] = h[['breaks']][1:22]
    plot(h, xlim=c(0, 21),
         col=color, border='white', cex.axis=1.5, cex.lab=1.5,
         xlab='Number of genes', ylab='', main='',
         yaxs='i', xaxs='i', lwd.axis=0, lwd.ticks=1, bty='l',
         xaxt='n', las=1, xpd=NA)
    axis(side=1, at=c(0, 5, 10, 15, 20.5),
         labels=c('0', '5', '10', '15', '>20'),
         cex.axis=1.5, lwd=0, lwd.ticks=1, gap.axis=0, xpd=NA)
    grid(nx=NA, ny=NULL)
    box(bty='l')

    hist(ld_data[ld_data[, 'assignment'] == replicon, 'pavp'],
         breaks=seq(0, 1, 0.05), xlim=c(0, 1),
         col=color, border='white', cex.axis=1.5, cex.lab=1.5,
         xlab='Proportion PAVs:total', ylab='', main='',
         yaxs='i', xaxs='i', lwd.axis=0, lwd.ticks=1, bty='l',
         las=1, gap.axis=0, xpd=NA)
    grid(nx=NA, ny=NULL)
    box(bty='l')
       
}

dev.off()
