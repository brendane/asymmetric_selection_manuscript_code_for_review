#!/usr/bin/env Rscript
#
# Compile relative fitnesses for the first experiment, the winter 2018
# experiments, the NPD experiment, and SpanFran1.
#
# No permutations -- that will be done separately.
#
# UPDATE 2020-10-26: added  in non-truncatula species experiment data.
#
# UPDATE 2021-01-09: Put all the SP1 results into one file.
#

library(magrittr)

PROJDIR = paste0(Sys.getenv('HOME'), '/project/select_reseq')
OUTDIR = file.path(PROJDIR, 'results/relative_fitness/compile/2020-05-23')

## Copy script to output directory
cargs = commandArgs(trailingOnly=FALSE)
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUTDIR, 'script_copies'), recursive=TRUE, showWarnings=FALSE)
file.copy(gsub('--file=', '', cargs[grepl('--file=', cargs)]),
          file.path(OUTDIR, 'script_copies'))



## For every experiment
## 1. Read fitnesses calculated previously
## 2. Alter strain names
## 3. (SpanFran1) split into medicae and meliloti
## 4. Write output table

## 101 strain community, first experiment
indir = file.path(PROJDIR, 'results/relative_fitness/160908_161228_reads/2018-08-16/real')
c101 = data.frame('strain'=character(101),
                  'C101_A17_Fit_mean'=numeric(101),
                  'C101_R108_Fit_mean'=numeric(101),
                  'C101_soil8wk_Fit_mean'=numeric(101),
                  'C101_TY72_Fit_mean'=numeric(101))
for(x in list(list('C101_A17_Fit_mean', 'bulk-HM101.tsv'),
              list('C101_R108_Fit_mean', 'bulk-HM340.tsv'),
              list('C101_soil8wk_Fit_mean', 'TY24-FS8W.tsv'),
              list('C101_TY72_Fit_mean', 'bulk-TY72.tsv'))) {
    coln = x[[1]]
    fname = x[[2]]
    fits = read.csv(file.path(indir, fname),
                            sep='\t', header=FALSE, as.is=TRUE)
    fits[, 1] = ifelse(fits[, 1] == '3011', fits[, 1], paste0('USDA', fits[, 1]))
    if(c101[1, 'strain'] == '') {
        c101[, 'strain'] = fits[, 1]
    } else {
        fits = fits[match(c101[, 'strain'], fits[, 1]), ]
    }
    c101[, coln] = fits[, 2]
}

## 101 strain experiment, multispecies
c101s = read.csv(file.path(PROJDIR, 'results/relative_fitness/101_strain_species/2020-10-26/real/phenotypes.tsv'),
                 sep='\t', header=TRUE, as.is=TRUE, check.names=FALSE)
names(c101s)[-1] = paste0("C101_", names(c101s)[-1])
c101 = merge(c101, c101s)

write.table(c101[order(c101[, 'strain']), ], file=file.path(OUTDIR, 'C101.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)



## 68 strain community, winter 2018 & NPD
w2018 = read.csv(file.path(PROJDIR, 'results/relative_fitness/winter_2018/2020-05-23/phenotypes.tsv'),
                 sep='\t', header=TRUE, as.is=TRUE, check.names=FALSE)
w2018[, 'strain'] = gsub('C$', 'c', w2018[, 'strain'])
w2018[, 'strain'] = ifelse(grepl('^[A-Z]', w2018[, 'strain']) | w2018[, 'strain'] == '3085',
                         w2018[, 'strain'],
                         paste0('USDA', w2018[, 'strain']))
scol = which(colnames(w2018) == 'strain')
names(w2018)[-scol] = paste0('C68_', names(w2018)[-scol])
npd = read.csv('results/relative_fitness/npd/2019-01-28/real/phenotypes.tsv',
               sep='\t', header=TRUE, as.is=TRUE, check.names=FALSE)
npd[, 'strain'] = ifelse(grepl('^[A-Z]', npd[, 'strain']) | npd[, 'strain'] == '3085',
                         npd[, 'strain'],
                         paste0('USDA', npd[, 'strain']))
scol = which(colnames(npd) == 'strain')
names(npd)[-scol] = paste0('NPD_', names(npd)[-scol])
c68 = merge(npd, w2018, by='strain')

write.table(c68[order(c68[, 'strain']), ], file=file.path(OUTDIR, 'C68.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)


## SpanFran1
mel_strains = scan(file.path(PROJDIR, 'data/strain/lists/spanfran1_meliloti.txt'),
                   what='character')
med_strains = scan(file.path(PROJDIR, 'data/strain/lists/spanfran1_medicae.txt'),
                   what='character')
sp1 = read.csv(file.path(PROJDIR, 'results/relative_fitness/spanfran1/2021-01-08_all_together/phenotypes.tsv'),
               sep='\t', header=TRUE, as.is=TRUE, check.names=FALSE)
sp1 = cbind(strain=sp1[, 'strain'],
            species=ifelse(sp1[,'strain'] %in% mel_strains, 'meliloti', 'medicae'),
            sp1[, -1])

write.table(sp1[order(sp1[, 'species'], sp1[, 'strain']), ],
            file=file.path(OUTDIR, 'SP1.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)


## SpanFran2 GWAS
sp2 = read.csv(file.path(PROJDIR, 'results/relative_fitness/spanfran2/2020-10-27/real/phenotypes.tsv'),
               sep='\t', header=TRUE, check.names=FALSE, as.is=TRUE)
scol = which(colnames(sp2) == 'strain')
names(sp2)[-scol] = paste0('SP2_', names(sp2)[-scol])
write.table(sp2, file=file.path(OUTDIR, 'SP2-gwas.tsv'), sep='\t',
            col.names=TRUE, row.names=FALSE, quote=FALSE)


## SpanFran NCR 86 strain
ncr = read.csv(file.path(PROJDIR, 'results/relative_fitness/ncr/2021-02-21/real/phenotypes.tsv'),
               sep='\t', header=TRUE, check.names=FALSE, as.is=TRUE)
scol = which(colnames(ncr) == 'strain')
names(ncr)[-scol] = paste0('NCR_', names(ncr)[-scol])
write.table(ncr, file=file.path(OUTDIR, 'NCR-C86.tsv'), sep='\t',
            col.names=TRUE, row.names=FALSE, quote=FALSE)
