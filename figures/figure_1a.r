#!/usr/bin/env Rscript

source('colors.r')

values = c('Shannon'=37,
           'Nod. number'=42,
           'Nod. area'=67,
           'Biomass'=49)

pdf('medicago_pve.2021-09-24.pdf', width=2.5, height=3.7)

barplot(values,
        xaxs='i', yaxs='i', cex.axis=1.5, cex.lab=1.5, main='',
        xaxt='n',
        ylab='PVE by genotype',
        xlab='', las=2,
        col=phenotype_colors[tolower(gsub('\\. ', '_', names(values)))],
        border=NA)
axis(side=1, at=1:length(values), cex.lab=2, cex.axis=1.5, labels=names(values),
     lwd.tick=0, lwd=0, gap.axis=0, las=2)

dev.off()
