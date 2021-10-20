#!/bin/bash
#
# Use FreeBayes to call SNPs and small indels on the meliloti strains
# collected by Mike Grillo using reads aligned to Rm41. Meant for
# comparison with JGI SNP calls.
#
# INPUT
#
# Reads trimmed with TrimGalore! / HT-Stream and aligned to Rm41.
# with bwa mem.
#
# UPDATE 26 October 2018: Switched to Rm41 assembly from Bielefeld
#    University to match JGI (this changed the output directory).
#
# SETTINGS
#
QUEUE="small"
WALLTIME="16:00:00"
NTHREADS=24
MEM="62GB"
#
REGIONSIZE=50000 # Size of regions to split the reference into for parallel operation
MAPQ=30          # Minimum mapping quality required for reads

SPECIES="MAG_meliloti"
REFSTRAIN="Rm41_bielefeld"
PROJDIR="${HOME}/project/jgi_sequencing"
INDIR="${PROJDIR}/results/align/${SPECIES}/${REFSTRAIN}/bwa/2018-10-11"
OUTDIR="${PROJDIR}/results/variant_calls/snps/${SPECIES}/${REFSTRAIN}/freebayes/2018-10-11"
REFERENCE="${PROJDIR}/data/reference/genome/${REFSTRAIN}/genome.fasta"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"
SCRATCHDIR="/scratch.global/${USER}/variant_calls_snps_${SPECIES}_${REFSTRAIN}_freebayes_2018-10-11"

vcf2tsv="${PROJDIR}/script/bin/vcf2tsv.py"
vcf2dgrp="${PROJDIR}/script/bin/vcf2dgrp_whole.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bcftools/1.6
    module load freebayes/20181002
    module load htslib/1.6
    module load parallel/20160822
    module load vcflib/20181002

    cd "$SCRATCHDIR"

    SECONDS=0

    fasta_generate_regions.py "${REFERENCE}.fai" "$REGIONSIZE" > \
        "ref.regions" \
        || { echo "generating list of regions failed"; exit 1; }

    rm -f "local_bams.txt"
    IFS=$'\\n' ls "$INDIR" | grep -v "work\\|array\\|log\\|script" > "input_files.txt" \
        || { echo "listing strains failed"; exit 1; }
    while read -r d; do
        local_name="$(pwd)/$(basename "$d").cram"
        cp "${INDIR}/${d}/alignment.cram" "$local_name" \
            || { echo "copying ${d} cram file failed"; exit 1; }
        cp "${INDIR}/${d}/alignment.cram.crai" "${local_name}.crai" \
            || { echo "copying ${d} crai file failed"; exit 1; }
        echo "$local_name" >> "local_bams.txt"
    done < "input_files.txt"

    freebayes-parallel "ref.regions" "$NTHREADS" \
        -f "$REFERENCE" -p 1 -L "local_bams.txt" \
        --min-mapping-quality "$MAPQ" > "tmp.vcf" \
            || { echo "calling variants failed"; exit 1; }

    # Compress and index the FreeBayes output -- needed for reading
    # with bcftools
    rm -f "tmp.vcf.gz"
    bgzip -i "tmp.vcf" \
        || { echo "gzipping tmp file failed"; exit 1; }
    tabix "tmp.vcf.gz" \
        || { echo "tabixing tmp file failed"; exit 1; }

    # Add variant IDs, change name, and compress and index
    bcftools annotate --set-id 'freebayes-snp-%CHROM-%POS-%FIRST_ALT' "tmp.vcf.gz" > \
        "variants.vcf" \
        || { echo "normalizing SNPs failed"; exit 1; }
    rm -f "variants.vcf.gz"
    bgzip -i "variants.vcf" \
        || { echo "gzipping variants failed"; exit 1; }
    tabix "variants.vcf.gz" \
        || { echo "tabixing variants failed"; exit 1; }

    # Statistics at different quality filter stringencies
    vcfstats "variants.vcf.gz" > "stats.all.txt" \
        || { echo "getting stats failed"; exit 1; }
    for q in $(seq 5 5 30); do
        vcffilter -f "QUAL > ${q}" "variants.vcf.gz" | \
            vcfstats > "stats.qual${q}.txt" \
            || { echo "getting stats for qual > ${q} failed"; exit 1; }
    done

    # Filter variants by quality score, then compress and index
    vcffilter -f "QUAL > 20" "variants.vcf.gz" > \
        "variants.qual20.vcf" \
        || { echo "filtering at QUAL >= 20 failed"; exit 1; }
    rm -f "variants.qual20.vcf.gz"
    bgzip -i "variants.qual20.vcf" \
        || { echo "gzipping filtered variants failed"; exit 1; }
    tabix "variants.qual20.vcf.gz" \
        || { echo "tabixing filtered variants failed"; exit 1; }

    # Make a tsv file for SNPs, reported by site.
    # Do no filtering on genotyping rate.
    # First have to turn MNPs into SNPs and might as well get rid
    # of indels at the same time.
    # Added bcftools norm because I was getting repeat rows -- mostly
    # in one spot around 3,600,000 on the chromosome with a weird MNP.
    vcfallelicprimitives "variants.qual20.vcf.gz" | \
        vcfnoindels | \
        bcftools norm -f "$REFERENCE" -D | \
        bcftools annotate --set-id 'freebayes-snp-%CHROM-%POS-%FIRST_ALT' > \
            "snps.qual20.primitives.vcf" \
        || { echo "making indel-free primitives file failed"; exit 1; }
    cp "$vcf2tsv" .
    ./$(basename "$vcf2tsv") --rs --het-miss \
        --output "snps.qual20.primitives.tsv" \
        --min-gt 0.0 "snps.qual20.primitives.vcf" \
        || { echo "making tsv file for filtered variants failed"; exit 1; }

    # Make an input file for HARP
    cp "$vcf2dgrp" .
    ./$(basename "$vcf2dgrp") --output "snps.qual20.primitives.dgrp" \
        --reference "$REFERENCE" \
        --pad 1000 --min-gt 0.0 "snps.qual20.primitives.vcf" \
        || { echo "making DGRP file for filtered variants failed"; exit 1; }

    # Bgzip
    rm -f "snps.qual20.primitives.vcf.gz"
    bgzip -i "snps.qual20.primitives.vcf" \
        || { echo "gzipping variants failed"; exit 1; }
    tabix "snps.qual20.primitives.vcf.gz" \
        || { echo "tabixing variants failed"; exit 1; }

    # Report stats for SNPs
    vcfstats "snps.qual20.primitives.vcf.gz" > \
        "stats.snps.qual20.primitives.txt" \
        || { echo "stats for SNPs failed"; exit 1; }

    # Copy from scratch directory
    cp "variants.vcf.gz"* "$OUTDIR" \
        || { echo "copying vcf and indexes failed"; exit 1; }
    cp "snps.qual20.primitives.vcf.gz"* "$OUTDIR" \
        || { echo "copying snps and vcf indexes failed"; exit 1; }
    cp "snps.qual20.primitives.tsv" "$OUTDIR" \
        || { echo "copying snps tsv failed"; exit 1; }
    cp "snps.qual20.primitives.dgrp" "$OUTDIR" \
        || { echo "copying snps DGRP failed"; exit 1; }
    cp "snps.qual20.primitives.dgrp.contigs" "$OUTDIR" \
        || { echo "copying DGRP contigs file failed"; exit 1; }
    cp "input_files.txt" "$OUTDIR" \
        || { echo "copying input list failed"; exit 1; }
    cp stats.* "$OUTDIR" \
        || { echo "copying stats files failed"; exit 1; }
    cp $(basename $vcf2tsv) "$OUTDIR" \
        || { echo "copying vcf2tsv script failed"; exit 1; }
    cp $(basename $vcf2dgrp) "$OUTDIR" \
        || { echo "copying vcf2dgrp script failed"; exit 1; }

#    rm * || { echo "deleting tmp files failed"; exit 1; }

    echo "RUN TIME = ${SECONDS} seconds ($(($SECONDS/60)) minutes)"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$SCRATCHDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-$(basename $0)"
    cp "$0" "$sfile"

    cd "$WORKDIR"
    qsub2 -A "tiffinp" -W group_list="tiffinp" -q "$QUEUE" \
        -l "walltime=${WALLTIME},nodes=1:ppn=${NTHREADS},mem=${MEM}" \
        "$sfile"
    echo "$WORKDIR"

fi
