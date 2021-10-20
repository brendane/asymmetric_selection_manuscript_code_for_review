#!/usr/bin/env Rscript
#
# Genotype means, not adjusted for greenhouse position, inoculated
# only.
#
# Include biomass and the last height measurement.
#

PROJDIR = '/home/tiffinp/epste051/project'
phenotype_file = file.path(PROJDIR, 'select_reseq/results/phenotypes/gwas_phenotypes/spanfran2/2021-06-29_nod_shannon/phenotypes.tsv')

hist_color = 'gray20'
points_color = rgb(0.4, 0.4, 0.4, 0.8)
source('colors.r')

phenotypes = read.csv(phenotype_file, sep='\t', header=TRUE, as.is=TRUE)


upper_panel = function(x, y, ...) {
    #points(x, y, col=points_color, ...)
    cr = cor(x, y)
    usr = par("usr")
    inset = c(0, 0)
    insetx = inset[1L] * (usr[2L] - usr[1L])
    insety = inset[2L] * (usr[4L] - usr[3L])
    w = abs(usr[2L] - usr[1L])
    #if(cr > 0) {
    #    xx = usr[1L] + insetx
    #    adj = c(0, 1)
    #} else {
    #    xx = usr[2L] - insetx
    #    adj = c(1, 1)
    #}
    #yy = topright = usr[4L] - insety
    xx = mean(par('usr')[c(1,2)])
    yy = mean(par('usr')[c(3,4)])
    colr = colorRamp(c('red', 'gray70', 'blue'))(cr / 2 + 0.5)
    adj = c(0.5, 0.5)
    cex = abs(cr) * 3 + 2
    text(xx, yy, round(cr, 2),
         col=rgb(colr[1]/255, colr[2]/255, colr[3]/255),
         cex=cex, adj=adj)
}


phenotypes_for_plotting = phenotypes[, c('host_PC1', 'host_PC2', 'shannon',
                                         'nod_number', 'nod_area')]
names(phenotypes_for_plotting) = c('PC1', 'PC2', 'Shannon diversity', 'Nodule number',
                                   'Nodule area')
cm = cor(phenotypes_for_plotting)

panel.hist <- function(x, ...) {
    ## Directly from example function
    phenotype = NULL
    for(ph in colnames(phenotypes_for_plotting)) {
        ph2 = gsub('pc', 'PC', gsub('nodule', 'nod',
                                    tolower(gsub(' ', '_', gsub('host_', '', gsub(' diversity', '', ph))))))
        if(all(phenotypes_for_plotting[, ph] == x)) {
            phenotype = ph2
            break
        }
    }
    colr = phenotype_colors[phenotype]
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE, breaks=7)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col=colr, border='white', ...)
}

pdf('plant_phenotype_correlations.2021-09-24.pdf', height=9, width=9)
pairs(phenotypes_for_plotting,
      diag.panel=panel.hist, upper.panel=NULL,
      lower.panel=upper_panel,
      cex.labels=2, cex.axis=1.5, cex.lab=1.5, las=1, bty='n')
dev.off()
