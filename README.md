# Files

- Table of orthologs for rhizobia:
    ensifer_variants/table_of_orthologs.tsv
 
- Rhizobia variant calls:
    Variant calls for the 88 strain community, used for association
    analyses: ensifer_variants/variant_calls_for_association.vcf.gz

- LD group assignments for rhizobia:
    ensifer_variants/ld_group_assignments.tsv

- Raw GEMMA output:
    To be included in the final Dryad package (too large for github).

- For Medicago variant calls please see www.medicagohapmap2.org.


# Code

These are the scripts used for processing the data. Note that they have
been uploaded "as-is", which means that they are designed to run in a
specific computing environment. As such, they are useful as a record of
how the data were processed, but would need significant modification
before they could be re-run.

## Medicago variant calls:

- Download, clean, and align Medicago reads:
    align_mt5.0_bwa_2020-05-08.sh

- Call SNPs and small indels by using FreeBayes:
    variant_calls_snps_freebayes_2020-07-03.sh

- Impute SNPs and small indels:
    impute_beagle_2020-08-17.sh

- Combine and filter variants for association mapping:
    combine_variants_spanfran2_medicago_2020-11-16_freebayes_rdvs.sh 
    (NOTE: This script refers to a set of read-depth variants, but these
     were filtered out of the final set of association analyses.)


## Choosing strains for the SpanFran community:

- Align reads to reference genome:
    align_MAG_meliloti_Rm41_bwa_2018-10-11.sh (used only for choosing strains)

- Variant calls (also only used for choosing strains)
    variant_calls_snps_MAG_meliloti_Rm41_freebayes_2018-10-11.sh

- Choose strains for the community:
    strain_choosing_code.r


## Handling individual rhizobial genome sequences:

- Clean reads:
    clean_reads_all_trimgalore_htstream_2018-09-25_dedup_phix_trim.sh

- Assemblies:
    Run by Alex at the Univ of IL

- Whole genome alignment and identification of orthologs:
    ortho_medmel_sibeliaz_2020-04-09_split.sh

- Alignment of orthologs using Pagan and identification of SNPs,
  also identification of medicae-meliloti cross-species orthologs:
    gene_aligns_from_sibeliaz_pagan2_2021-04-04.sh

- Replicon assignments for genes:
    replicon_assignments_2021-04-05.sh

- Combine and filter variants for association mapping and LD grouping:
    combine_variants_meliloti_2021-04-06_denovo.sh
    Note to self: This script does not use the "filtered" output from
    the ortholog step, but that is okay--the differences from the
    filtered output do not affect this step for the SpanFran dataset.

- LD groups:
    ld_pangenome_r2_groups_2021-04-06_denovo.sh


## Strain frequencies:

- Variants aligned to reference genome:
    clean_reads_all_trimgalore_htstream_2018-09-25_dedup_phix_trim.sh
        (listed under rhizobia variant calls)
    clean_reads_illinois_trimgalore_htstream_2020-04-15_dedup_phix_trim.sh
    align_meliloti_USDA1106_bwa_2018-10-11.sh
    align_meliloti_USDA1106_bwa_2020-04-16.sh
    variant_calls_snps_meliloti_USDA1106_freebayes_2020-04-17_with_33_extra.sh

- Clean reads:
    trim_reads_spanfran2_trim_galore_2020-10-23_gwas.sh

- Align to reference:
    align_spanfran2_bwa_2020-10-24_gwas.sh

- Estimate strain frequencies from core SNPs:
    strain_frequency_spanfran2_harp_2020-10-24_core.sh

- Relative fitness calculation:
    relative_fitness_spanfran2_2020-10-27.r
    relative_fitness_compile_2020-05-23.r


## Association analyses:

- PCA phenotypes:
    phenotypes_gwas_phenotypes_spanfran2_2020-11-13_pca.r 

- Shannon's index and nodule phenotypes:
    adjusted means from Peter Tiffin
    phenotypes_gwas_phenotypes_spanfran2_2021-06-29_nod_shannon.r

- GEMMA analysis:
    - Rhizobia PCs:
        genome_scan_spanfran2_gwas_2021-06-29_irnn_by_replicon.sh
    - Old rhizobia data:
        genome_scan_spanfran2_gwas_2021-07-13_old_data.sh
    - Medicago PCs:
        genome_scan_spanfran2_gwas_host_2020-11-18_lmm.sh
    - Medicago Shannon & Nodule phenotypes:
        genome_scan_spanfran2_gwas_host_2021-06-29_nod_traits_lmm.sh

- Host-allele-specific fitness:
    phenotypes_rhizobia_phenotypes_projected_frequency_2021-12-02.sh
    genome_scan_spanfran2_allele_specific_w_2021-04-07.sh
    genome_scan_spanfran2_allele_specific_w_2021-04-07_randomized.sh

    Note: final_freqs.tsv file is dated to June 7, while major results files
    are dated May 27. However, this appears to be just when the file was
    copied from the scratch directory, because there are no script runs
    anytime around June 7.


## Population genetics analysis:

- Ensifer popgen stats from libsequence: popgen_summaries_2021-04-05_libsequence.sh
