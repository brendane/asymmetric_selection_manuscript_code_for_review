#!/usr/bin/env Rscript
#
# Manhattan plots for SNPs in Medicago. Both effect size and p-values.
#
# Uses a strategy of vector graphics for axes and labels and raster graphics
# for the points, because there are too many points to comfortably display
# in vector format.
#
# LMM GEMMA runs shown.
#

WINDOW_FRACTION = 0.20

fai_file = '/home/tiffinp/epste051/project/mtr_variant_calls/data/assemblies/Mtr_genomes/MtrunA17r5.0-ANR/genome.fasta.fai'
data_file = '/home/tiffinp/epste051/project/select_reseq/results/summary_analyses_and_data_collections/analyses/spanfran2_gwas/gwas_analysis_2021-07-30/medicago_results_spreadsheets/medicago_data_compiled_by_window.tsv'
candidates_file = '/home/tiffinp/epste051/project/select_reseq/results/summary_analyses_and_data_collections/analyses/spanfran2_gwas/gwas_analysis_2021-07-30/candidate_lists/medicago_candidate_windows.tsv'
phenotypes = c('PC1'='host_score_PC1', 'PC2'='host_score_PC2', 'Shannon\'s Diversity'='shannon',
               'Nodule Area'='nod_area', 'Nodule Number'='nod_number')
chromosomes = c('CM010648.1'='Chr1',
                'CM010649.1'='Chr2',
                'CM010650.1'='Chr3',
                'CM010651.1'='Chr4',
                'CM010652.1'='Chr5',
                'CM010653.1'='Chr6',
                'CM010654.1'='Chr7',
                'CM010655.1'='Chr8')
pad_prop = 0.02

candidates = read.csv(candidates_file, sep='\t', header=TRUE, as.is=TRUE)
window_data = read.csv(data_file, sep='\t', comment.char='#', header=TRUE, as.is=TRUE)

fai = read.table(fai_file, sep='\t', as.is=TRUE, header=FALSE)
chr_lengths = fai[, 2]
pad_bp = sum(chr_lengths) * pad_prop
starts = structure(head(c(0, cumsum(chr_lengths + pad_bp)),n=-1),
                   names=fai[, 1])

do_manhattan_plot = function(window_data, fname, phenotype, candidate_windows) {
    w = 7.5
    h = 3

    ps = window_data[['start']] + starts[ window_data[['chromosome']] ]
    wins = paste(window_data[['chromosome']], window_data[['start']],
                 window_data[['end']], sep='-')

    ## Plot the axes, labels, etc.
    #pdf(fname, width=w, height=h)
    #par('mgp'=c(1.5, 0.5, 0))

    mx = max(ps, na.rm=TRUE)
    mn = min(ps, na.rm=TRUE)
    y = window_data[, paste0(phenotypes[phenotype],
                             '_max_effect_size_lmm')]
    i = rank(-y) <= length(y) * WINDOW_FRACTION
    y = y[i]
    ps = ps[i]
    wins = wins[i]

    plot(ps, y,
         main=phenotype, ylab='Effect size',
         xlab='Genomic position',
         xaxt='n', bty='L', xaxs='i', yaxs='i',
         cex=0.5, xlim=c(mn, mx),
         col=ifelse(wins %in% candidate_windows, 'red', rgb(0.2, 0.2, 0.2, 0.4)),
         xpd=NA, las=1)
    cnd_points = cbind(ps[wins %in% candidate_windows],
                       y[wins %in% candidate_windows])

    ## Axis
    axis(side=1, at=mapply(function(x, y) mean(c(x,y)),
                           fai[1:8, 2]+starts[1:8],
                           starts[1:8]),
         tick=FALSE,
         labels=chromosomes[names(starts[1:8])],
         gap.axis=0.1)

    #dev.off()
}

pdf('medicago_manhattan_effect_size_all.pdf', width=7.5, height=12)
par('mgp'=c(0.5, 0.5, 0), 'mfrow'=c(5, 1))
for(phenotype in names(phenotypes)) {

    candidate_windows = candidates[candidates[, 1] == gsub('host_score_', '', phenotypes[phenotype]) &
                                   candidates[, 2] == 'effect_size' &
                                   candidates[, 3] == 0.01 &
                                   candidates[, 4] == 'candidate',
                                   5]

    do_manhattan_plot(window_data,
                      paste0('medicago_manhattan_effect_size_', phenotype, '.pdf'),
                      phenotype,
                      candidate_windows)

}
