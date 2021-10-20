#!/usr/bin/env Rscript
#
# Run PCA with centering but not scaling on median fitness matrices.
# Report both the scores (the "x" component of the prcomp output) and
# the loadings ("rotation" component).
#
# Limit to just the truncatula strains.
#
# UPDATE 2021-07-08: Make a table of individual shannon index values.
#
# UPDATE 2021-08-13: Added Shannon's equitability.
#
# NOTE: To use the spatialEco package, need to load proj, geos, udunits,
# and gdal modules.
#

library(vegan)
library(spatialEco)

PROJDIR = paste0(Sys.getenv('HOME'), '/project/select_reseq')
OUTDIR = file.path(PROJDIR, 'results/phenotypes/gwas_phenotypes/spanfran2/2021-06-29_nod_shannon')

infile = file.path(PROJDIR, 'results/relative_fitness/compile/2020-05-23/SP2-gwas.tsv')
rep_freq_file = file.path(PROJDIR, 'results/strain_frequency/spanfran2/harp/2020-10-24_core/freq.tsv')
ph_file = file.path(PROJDIR, 'results/phenotypes/mtr_phenotypes/spanfran2/from_peter_2020-08-11/GWAS2_pheno_aug2020_plt.csv')
means_file = file.path(PROJDIR, 'results/phenotypes/mtr_phenotypes/spanfran2/from_peter_2021-06-28/means.tsv')

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
med_fit = as.matrix(X[, grepl('_Freq_med', colnames(X))])
colnames(med_fit) = gsub('SP2_', '', gsub('_Freq_med', '', colnames(med_fit)))
rownames(med_fit) = X[, 'strain']
med_fit = t(med_fit)[truncatula, ]

## Shannon's diversity
shannon = cbind('strain'=rownames(med_fit),
                'shannon'=diversity(med_fit, MARGIN=1, index='shannon'))
## Add a little bit to the frequencies to avoid zeroes (there is a log
## operation). Also, treat as counts because the median frequencies don't
## always sum to 1.
shannon_e = cbind('strain'=rownames(med_fit),
                  'evenness'=shannons(med_fit+1E-5, margin='row', ens=FALSE, counts=TRUE)[, 'evenness'])


## Nod morphology
morph = read.csv(means_file, sep='\t', header=TRUE, as.is=TRUE)
names(morph)[names(morph) == 'geno'] = 'strain'

## combine
output = merge(merge(shannon, morph, all.x=TRUE), shannon_e, all=TRUE)


## Write output
write.table(output,
            file=file.path(OUTDIR, 'phenotypes.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

rep_freq = read.csv(rep_freq_file, sep='\t', header=TRUE, as.is=TRUE)
names(rep_freq)[names(rep_freq) == 'pool'] = 'plant'
geno_table = read.csv(ph_file, sep=',', header=TRUE, as.is=TRUE)[, c('plant', 'geno')]
shannon_by_rep = cbind('plant'=rep_freq[, 'plant'],
                       'shannon'=diversity(rep_freq[, -1], MARGIN=1, index='shannon'))
shannon_e_by_rep = cbind('plant'=rep_freq[, 'plant'],
                         'shannon_e'=shannons(rep_freq[, -1]+0.00001, margin='row',
                                              ens=FALSE, counts=FALSE)[, 'evenness'])
shannon_by_rep = merge(shannon_by_rep, shannon_e_by_rep, all=TRUE)
write.table(merge(shannon_by_rep, geno_table),
            file=file.path(OUTDIR, 'shannon_replicates.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
