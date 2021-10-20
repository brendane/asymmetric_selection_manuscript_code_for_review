#!/bin/bash
#
# Use FreeBayes to call SNPs and small indels on all the meliloti strains
# using reads aligned to USDA1106.
#
# Note: It might be desirable to set --no-population-priors in the future,
# however, I did not do that in this script. I wanted to keep the settings
# consistent with the medicae analysis.
#
# INPUT
#
# Reads trimmed with TrimGalore! / HT-Stream and aligned to USDA 1106
# with bwa mem.
#
# SETTINGS
#
QUEUE=amdsmall
WALLTIME=48:00:00
MEM=8GB
NTHREADS=1
#
REGIONSIZE=50000  # Size of regions to split the reference into for parallel operation
MAPQ=30           # Minimum mapping quality required for reads
N_ALLELES=4       # No more than this many alleles per variant (to save memory--only for gvcf)

SPECIES=meliloti
REFSTRAIN=USDA1106
PROJDIR="${HOME}/project/jgi_sequencing"
INDIR1="${PROJDIR}/results/align/${SPECIES}/${REFSTRAIN}/bwa/2018-10-11"
INDIR2="${PROJDIR}/results/align/${SPECIES}/${REFSTRAIN}/bwa/2020-04-16"
OUTDIR="${PROJDIR}/results/variant_calls/snps/${SPECIES}/${REFSTRAIN}/freebayes/2020-04-17_with_33_extra"
REFERENCE="${PROJDIR}/data/reference/genome/${REFSTRAIN}/genome.fasta"
OLDPYTHONPATH="$HOME/usr/lib/python2.7/site-packages.old"
SPANFRAN="${PROJDIR}/../select_reseq/data/strain/lists/meliloti_sibeliaz_with_freebayes_2020-04-09.txt"

WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYJOBDIR="${OUTDIR}/arraydata"
SCRATCHDIR="/scratch.global/${USER}/variant_calls_snps_${SPECIES}_${REFSTRAIN}_freebayes_2020-04-17"
TEMPVCFDIR="${SCRATCHDIR}/tmpvcf"

vcf2tsv="${PROJDIR}/script/bin/vcf2tsv.py"
vcf2dgrp="${PROJDIR}/script/bin/vcf2dgrp_whole.py"

if [[ "$PBS_ENVIRONMENT" == PBS_BATCH ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bcftools/1.9
    module load freebayes/20190623
    module load htslib/1.9
    module load parallel/20190122
    module load vcflib/1.0.1
    module unload python

    SECONDS=0

    TASK="$(cut -f 1 "${ARRAYJOBDIR}/${SLURM_ARRAY_TASK_ID}")"
    REGION="$(cut -f 2 "${ARRAYJOBDIR}/${SLURM_ARRAY_TASK_ID}")"

    cd "$SCRATCHDIR"

    case "$TASK" in

        "prep")
            rm -f "local_bams.txt"
            IFS=$'\n' ls "$INDIR1/" | awk '{ print "'"$INDIR1"'/" $1 };' | \
                grep -v "work\|array\|log\|script" > "input_files.txt" \
                || { echo "listing strains failed"; exit 1; }
            IFS=$'\n' ls "$INDIR2/" | awk '{ print "'"$INDIR2"'/" $1 };' | \
                grep -v "work\|array\|log\|script" >> "input_files.txt" \
                || { echo "listing strains failed"; exit 1; }
            while read -r d; do
                local_name="$(pwd)/$(basename "$d").cram"
                cp "${d}/alignment.cram" "$local_name" \
                    || { echo "copying ${d} cram file failed"; exit 1; }
                cp "${d}/alignment.cram.crai" "${local_name}.crai" \
                    || { echo "copying ${d} crai file failed"; exit 1; }
                echo "$local_name" >> "local_bams.txt"
            done < "input_files.txt"
            ;;

        "run")
            mkdir -p "$TEMPVCFDIR"
            rm -f "${TEMPVCFDIR}/tmp.${REGION}.vcf"
#            freebayes \
#                --region "$REGION" \
#                -f "$REFERENCE" \
#                -p 1 \
#                -L "local_bams.txt" \
#                --min-mapping-quality "$MAPQ" > \
#                "${TEMPVCFDIR}/tmp.${REGION}.vcf" \
#                || { echo "calling variants for ${REGION} failed"; exit 1; }
            freebayes \
                --region "$REGION" \
                --gvcf \
                -f "$REFERENCE" \
                -p 1 \
                --use-best-n-alleles "$N_ALLELES" \
                -L "local_bams.txt" \
                --min-mapping-quality "$MAPQ" > \
                "${TEMPVCFDIR}/tmp.${REGION}.gvcf" \
                || { echo "calling variants GVCF for ${REGION} failed"; exit 1; }
            ;;

        "combine")
            ## Concatenate all the tmp files using the same strategy
            ## employed by freebayes-parallel
#            cat $(ls "$TEMPVCFDIR/"*.gvcf | sort -t: -k 1b,1 -k 2n,2) \
#                | vcffirstheader \
#                | vcfstreamsort -w 1000 \
#                | vcfuniq \
#                > "merged.gvcf" \
#                || { echo "combining vcfs failed"; exit 1; }
#            bgzip -i -f merged.gvcf \
#                || { echo "gzipping gvcf failed"; exit 1; }
            tabix -p vcf merged.gvcf.gz \
                || { echo "tabixing gvcf failed"; exit 1; }

#            cat $(ls "$TEMPVCFDIR/"* | sort -t: -k 1b,1 -k 2n,2) | \
#                vcffirstheader | \
#                vcfstreamsort -w 1000 | \
#                vcfuniq > \
#                "tmp.vcf" \
#                || { echo "combining vcfs failed"; exit 1; }
#
#            ## Prepare for use with bcftools
#            rm -f "tmp.vcf.gz"
#            bgzip -i "tmp.vcf" \
#                || { echo "gzipping tmp file failed"; exit 1; }
#            tabix "tmp.vcf.gz" \
#                || { echo "tabixing tmp file failed"; exit 1; }
#
#            # Add variant IDs, change name, and compress and index
#            bcftools annotate --set-id 'freebayes-snp-%CHROM-%POS-%FIRST_ALT' "tmp.vcf.gz" > \
#                "variants.vcf" \
#                || { echo "normalizing SNPs failed"; exit 1; }
#            rm -f "variants.vcf.gz"
#            bgzip -i "variants.vcf" \
#                || { echo "gzipping variants failed"; exit 1; }
#            tabix "variants.vcf.gz" \
#                || { echo "tabixing variants failed"; exit 1; }
#
#            # Statistics at different quality filter stringencies
#            vcfstats "variants.vcf.gz" > "stats.all.txt" \
#                || { echo "getting stats failed"; exit 1; }
#            for q in $(seq 5 5 30); do
#                vcffilter -f "QUAL > ${q}" "variants.vcf.gz" | \
#                    vcfstats > "stats.qual${q}.txt" \
#                    || { echo "getting stats for qual > ${q} failed"; exit 1; }
#            done
#
#            # Filter variants by quality score, then compress and index
#            vcffilter -f "QUAL > 20" "variants.vcf.gz" > \
#                "variants.qual20.vcf" \
#                || { echo "filtering at QUAL >= 20 failed"; exit 1; }
#            rm -f "variants.qual20.vcf.gz"
#            bgzip -i "variants.qual20.vcf" \
#                || { echo "gzipping filtered variants failed"; exit 1; }
#            tabix "variants.qual20.vcf.gz" \
#                || { echo "tabixing filtered variants failed"; exit 1; }
#
#            # Make a tsv file for SNPs, reported by site.
#            # Do no filtering on genotyping rate.
#            # First have to turn MNPs into SNPs and might as well get rid
#            # of indels at the same time. norm -D gets rid of some repeat
#            # rows that are created by this process.
#            vcfallelicprimitives "variants.qual20.vcf.gz" | \
#                vcfnoindels | \
#                bcftools norm -f "$REFERENCE" -d "any" | \
#                bcftools annotate --set-id 'freebayes-snp-%CHROM-%POS-%FIRST_ALT' > \
#                "snps.qual20.primitives.vcf" \
#                || { echo "making indel-free primitives file failed"; exit 1; }
#            cp "$vcf2tsv" .
#            PYTHONPATH="$OLDPYTHONPATH" ./$(basename "$vcf2tsv") --rs --het-miss \
#                --output "snps.qual20.primitives.tsv" \
#                --min-gt 0.0 "snps.qual20.primitives.vcf" \
#                || { echo "making tsv file for filtered variants failed"; exit 1; }
#
#            # Make an input file for HARP
#            cp "$vcf2dgrp" .
#            PYTHONPATH="$OLDPYTHONPATH" ./$(basename "$vcf2dgrp") --output "snps.qual20.primitives.dgrp" \
#                --reference "$REFERENCE" \
#                --pad 1000 --min-gt 0.0 "snps.qual20.primitives.vcf" \
#                || { echo "making DGRP file for filtered variants failed"; exit 1; }
#
#            # Bgzip
#            rm -f "snps.qual20.primitives.vcf.gz"
#            bgzip -i "snps.qual20.primitives.vcf" \
#                || { echo "gzipping variants failed"; exit 1; }
#            tabix "snps.qual20.primitives.vcf.gz" \
#                || { echo "tabixing variants failed"; exit 1; }
#
#            # Report stats for SNPs
#            vcfstats "snps.qual20.primitives.vcf.gz" > \
#                "stats.snps.qual20.primitives.txt" \
#                || { echo "stats for SNPs failed"; exit 1; }
#
#            # Copy from scratch directory
#            cp "variants.vcf.gz"* "$OUTDIR" \
#                || { echo "copying vcf and indexes failed"; exit 1; }
#            cp "snps.qual20.primitives.vcf.gz"* "$OUTDIR" \
#                || { echo "copying snps and vcf indexes failed"; exit 1; }
#            cp "snps.qual20.primitives.tsv" "$OUTDIR" \
#                || { echo "copying snps tsv failed"; exit 1; }
#            cp "snps.qual20.primitives.dgrp" "$OUTDIR" \
#                || { echo "copying snps DGRP failed"; exit 1; }
#            cp "snps.qual20.primitives.dgrp.contigs" "$OUTDIR" \
#                || { echo "copying DGRP contigs file failed"; exit 1; }
#            cp "input_files.txt" "$OUTDIR" \
#                || { echo "copying input list failed"; exit 1; }
#            cp stats.* "$OUTDIR" \
#                || { echo "copying stats files failed"; exit 1; }
#            cp $(basename $vcf2tsv) "$OUTDIR" \
#                || { echo "copying vcf2tsv script failed"; exit 1; }
#            cp $(basename $vcf2dgrp) "$OUTDIR" \
#                || { echo "copying vcf2dgrp script failed"; exit 1; }
            cp merged.gvcf* "$OUTDIR" \
                || { echo "copying gvcf failed"; exit 1; }

            ;;

        spanfran)
            cd "$OUTDIR"
            ## Make a file with just SpanFran segregating sites to
            ## compare with gene alignments

            sed 's/__.\+//g' "$SPANFRAN" \
                | tr '[:lower:]' '[:upper:]' \
                | sed 's/MEL/mel/g' | sed 's/MED/med/g' \
                | sed 's/_/-/g' \
                | grep "^MAG" \
                > spanfran.txt \
                || { echo "making list of strains failed"; exit 1; }

            infile=snps.qual20.primitives.vcf.gz
            outfile=spanfran_snps.qual20.primitives.vcf
            bcftools view -Ov -S spanfran.txt "$infile" \
                | bcftools view -Ov --min-ac 1:nonmajor \
                > "$outfile" \
                || { echo "filtering failed"; exit 1; }
            bgzip -i -f "$outfile" \
                || { echo "bgzipping failed"; exit 1; }
            tabix -p vcf "${outfile}.gz" \
                || { echo "tabixing failed"; exit 1; }
            ;;

    esac

    rm "${ARRAYJOBDIR}/${SLURM_ARRAY_TASK_ID}"

    echo "RUN TIME = ${SECONDS} seconds ($(($SECONDS/60)) minutes)"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYJOBDIR"
    mkdir -p "$SCRATCHDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-$(basename $0)"
    cp "$0" "$sfile"
    module load freebayes/20190623

    i=0
    case "$1" in
        "prep")
            WALLTIME=00:15:00
            MEM=2GB
            echo "$1"$'\t'"_" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;
        "run")
            rm -rf "$TEMPVCFDIR"
            for region in $(fasta_generate_regions.py "${REFERENCE}.fai" "$REGIONSIZE"); do
                echo "$1"$'\t'"$region" > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            ;;
        "combine")
            WALLTIME=20:00:00
            echo "$1"$'\t'"_" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;
        "spanfran")
            WALLTIME=02:00:00
            echo "$1"$'\t'"_" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;
    esac

    if [[ "$i" == 1 ]]; then
        array=0
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    echo "$WORKDIR"
    sbatch --account=tiffinp --export=PBS_ENVIRONMENT=PBS_BATCH \
        --mem=$MEM --partition=$QUEUE --nodes=1 --ntasks-per-node=$NTHREADS \
        --time=$WALLTIME --array="$array" "$sfile"
    echo "Submitted ${i} jobs"
fi
