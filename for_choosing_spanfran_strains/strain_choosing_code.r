#!/usr/bin/env Rscript

## This code copied from notebook entry 2019-04-02.

threshold = 5000

## Read SNPs
meliloti_snps_file = '../../jgi_sequencing/results/variant_calls/snps/MAG_meliloti/Rm41_bielefeld/freebayes/2018-10-11/snps.qual20.primitives.tsv'
meliloti_snps = as.matrix(fread(meliloti_snps_file)[, -(1:3)])
meliloti_phenotyped = paste0('MAG',
                             scan('../data/strain/lists/mag_phenotyped_2019-03-27.txt',
                                  what='character'))

## Work with only SNPs that are genotyped in every strain, because
## that's what we'll use to estimate strain frequencies.
meliloti_mag_snps = meliloti_snps[, grepl('^MAG', colnames(meliloti_snps))]
meliloti_mag_snps_gt = apply(meliloti_mag_snps, 1, function(x) sum(!is.na(x)))
meliloti_mag_snps_allgt = meliloti_mag_snps[meliloti_mag_snps_gt == ncol(meliloti_mag_snps), ]

## Count the number of differences
meliloti_ndiff = matrix(nrow=ncol(meliloti_mag_snps),
                        ncol=ncol(meliloti_mag_snps),
                        dimnames=list(colnames(meliloti_mag_snps),
                                      colnames(meliloti_mag_snps)),
                        data=NaN)
for(i in 1:(ncol(meliloti_ndiff)-1)) {
    s1 = colnames(meliloti_ndiff)[i]
    for(j in (i+1):ncol(meliloti_ndiff)) {
        s2 = colnames(meliloti_ndiff)[j]
        nd = sum(meliloti_mag_snps_allgt[, s1] != meliloti_mag_snps_allgt[, s2])
        meliloti_ndiff[s1, s2] = nd
        meliloti_ndiff[s2, s1] = nd
    }
}

m = meliloti_ndiff[rownames(meliloti_ndiff) %in% meliloti_phenotyped,
                   colnames(meliloti_ndiff) %in% meliloti_phenotyped]

## Assemble a community
n_strains = structure(numeric(30), names=(1:30)*200)
strains1000 = NULL
strains5000 = NULL
for(threshold in as.numeric(names(n_strains))) {
    strain_list = colnames(m)
    exclude_strains = logical(length(strain_list))
    names(exclude_strains) = strain_list
    number_above_threshold = apply(m, 1,
                                   function(x) sum(x > threshold, na.rm=TRUE))
    mean_simm = colMeans(m)
    stopifnot(all(strain_list == colnames(m)))
    for(s1 in strain_list[order(number_above_threshold,
                                mean_simm, decreasing=TRUE)]) {
        for(s2 in strain_list[order(number_above_threshold,
                                    mean_simm, decreasing=TRUE)]) {
            if(s1 == s2) {
                next
            }
            if(exclude_strains[s1] || exclude_strains[s2]) {
                next
            }
            if(m[s1, s2] < threshold) {
                s1_nat = number_above_threshold[s1]
                s2_nat = number_above_threshold[s1]
                if(s1_nat < s2_nat) {
                    exclude_strains[s1] = TRUE
                } else if(s2_nat < s1_nat) {
                    exclude_strains[s2] = TRUE
                } else {
                    if(s1 %in% meliloti_phenotyped && !(s2 %in% meliloti_phenotyped)) {
                        exclude_strains[s1] = TRUE
                    } else if(s2 %in% meliloti_phenotyped && !(s1 %in% meliloti_phenotyped)) {
                        exclude_strains[s2] = TRUE
                    } else {
                        exclude_strains[s1] = TRUE
                        warning('picked a strain arbitrarily')
                    }
                }
            }
        }
    }
    simple_strains = strain_list[!exclude_strains]
    n_strains[as.character(threshold)] = length(simple_strains)
    if(threshold == 1000) {
        strains1000 = simple_strains
    } else if(threshold == 5000) {
        strains5000 = simple_strains
    }
}

cat(strains1000, sep='\n', file='table/meliloti_mag_strains_chosen_1000_2019-04-02.txt')
cat(strains5000, sep='\n', file='table/meliloti_mag_strains_chosen_5000_2019-04-02.txt')

## Note: including strains that have /not/ been phenotyped adds two strains
## to both the 1,000 difference and 5,000 difference communities.
