#!/usr/bin/env Rscript
#
# Compile relative fitness values and other related traits for the SpanFran2
# experiment samples.
#

MIN_FIT_VALUE = -8
N_RANDOM = 1000

projdir = file.path(Sys.getenv('HOME'), 'project', 'select_reseq')
trt_file = file.path(projdir, 'results/phenotypes/mtr_phenotypes/spanfran2/from_peter_2020-08-11/GWAS2_pheno_aug2020_plt.csv')
outdir = file.path(projdir, 'results', 'relative_fitness', 'spanfran2',
                   '2020-10-27')
indirs = c(file.path(projdir, 'results', 'strain_frequency', 'spanfran2',
                     'harp', '2020-10-24_core'))
#           file.path(projdir, 'results', 'strain_frequency', 'spanfran2',
#                     'harp', '2020-10-27_ncr'))

cargs = commandArgs(trailingOnly=FALSE)
script_name = gsub('--file=', '', cargs[grepl('--file=', cargs)])

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir, 'random'), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir, 'real'), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir, 'script_copies'), showWarnings=FALSE)
file.copy(script_name, file.path(outdir, 'script_copies',
                                 paste0(format(Sys.time(), '%Y-%m-%d-%H%M'),
                                        '-', basename(script_name))))

trts = read.csv(trt_file, sep=',', as.is=TRUE, check.names=FALSE,
                header=TRUE)

fs = vector('list', length=length(indirs))
for(i in seq_along(indirs)) {
    f = read.csv(file.path(indirs[i], 'freq.tsv'),
                 sep='\t', header=TRUE, as.is=TRUE, check.names=FALSE)
    if(i > 1) {
        if(!all(colnames(f) %in% colnames(fs[[1]]))) {
            stop('Incompatible files')
        }
        f = f[, colnames(fs[[1]])]
    }

    fs[[i]] = f
}

all_freqs = do.call(rbind, fs)
freqs = as.matrix(all_freqs[, -1])

pool = f[, 'pool']
host = trts[match(f[, 'pool'], trts[, 'plant']), 'geno']
host[grepl('_I_', pool)] = 'initial'

mean_init_freq = colMeans(freqs[host == 'initial', ])
rep_fits = apply(freqs[host != 'initial', ], 1, function(x) log2(x / mean_init_freq))
rep_fits[is.infinite(rep_fits)] = MIN_FIT_VALUE
colnames(rep_fits) = all_freqs[host != 'initial', 'pool']

host_names = unique(host[host != 'initial'])
med_fits = matrix(nrow=ncol(freqs), ncol=length(host_names),
                  dimnames=list(rownames(rep_fits), host_names))
mean_fits = matrix(nrow=ncol(freqs), ncol=length(host_names),
                   dimnames=list(rownames(rep_fits), host_names))
med_freq = matrix(nrow=ncol(freqs), ncol=length(host_names),
                  dimnames=list(rownames(rep_fits), host_names))
for(h in host_names) {
    mf = apply(freqs[host == h, ], 2, function(x) median(x))
    med_freq[rownames(rep_fits), h] = mf[rownames(rep_fits)]
    med_fits[rownames(rep_fits), h] =
    apply(rep_fits[, host[host != 'initial'] == h], 1, median)
    mean_fits[rownames(rep_fits), h] =
    apply(rep_fits[, host[host != 'initial'] == h], 1, mean)

    med_rands = matrix(nrow=nrow(med_fits), ncol=N_RANDOM,
                       dimnames=list(rownames(med_fits)))
    for(i in 1:N_RANDOM) {
        med_rands[rownames(med_fits), i] = sample(med_fits[, h], nrow(med_rands), replace=FALSE)
    }
    write.table(cbind('strain'=rownames(med_rands), med_rands),
                file=file.path(outdir, 'random', paste0(h, '_Fit_med.tsv')),
                sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
}

colnames(med_fits) = paste0(colnames(med_fits), '_Fit_med')
colnames(mean_fits) = paste0(colnames(mean_fits), '_Fit_mean')
colnames(med_freq) = paste0(colnames(med_freq), '_Freq_med')
stopifnot(all(rownames(med_freq) == rownames(mean_fits)) && all(rownames(med_fits) == rownames(mean_fits)))
output = cbind(med_fits, mean_fits, med_freq)

write.table(cbind('strain'=rownames(output), output),
            file=file.path(outdir, 'real', 'phenotypes.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

write.table(cbind('pot'=colnames(rep_fits),
                  'host'=trts[match(colnames(rep_fits), as.character(trts[, 'plant'])), 'geno'],
                  t(rep_fits)),
            file=file.path(outdir, 'real', 'fitness_by_replicate.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
