#!/bin/bash
#
# Use HARP (Kessner et al. 2013) to estimate strain frequency for the
# communities from the winter 2020 SpanFran2 GWAS experiment, using
# reads aligned to USDA1106.
#
# Require that sites be genotyped in all strains.
#
# UPDATE 2021-05-26: Changed the variant input file because the previous
# one was in the wrong order.
#
# SETTINGS
#
WALLTIME="02:30:00"
QUEUE=amdsmall
MEM="16GB"
NTHREADS=4

PROJDIR="${HOME}/project/select_reseq"
OUTDIR="${PROJDIR}/results/strain_frequency/spanfran2/harp/2020-10-24_core"
BAMDIR="${PROJDIR}/results/align/spanfran2/bwa/2020-10-24_gwas"
REFERENCE="${PROJDIR}/data/reference/USDA1106/USDA1106.final.fasta"
STRAIN_FILE="${PROJDIR}/data/strain/lists/spanfran2_meliloti.txt"
SNPCALLS="${PROJDIR}/../jgi_sequencing/results/variant_calls/snps/meliloti/USDA1106/freebayes/2020-04-17_with_33_extra/snps.qual20.primitives.dgrp"

SCRATCHDIR="/scratch.global/${USER}/strain_frequency_spanfran2"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
ARRAYDIR="${OUTDIR}/arrayjobdata"
SCRIPTDIR="${OUTDIR}/script_copies"

concat_fasta="${PROJDIR}/script/bin/concatenate_reference_fasta.py"
concat_bam="${PROJDIR}/script/bin/concatenate_reference_bam.py"
filter_dgrp="${PROJDIR}/script/bin/filter_dgrp.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load harp/130306
    module load samtools/1.9

    set -euo pipefail

    TASK="$(cut -f 1 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"
    POOL="$(cut -f 2 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"

    cd "$OUTDIR"

    case "$TASK" in

        "prep")

            cd "$SCRATCHDIR"

            cp "${STRAIN_FILE}" "keep.txt" \
                || { echo "copying strains file failed"; exit 1; }

            cut -f 1  "${SNPCALLS}.contigs" \
                | sed 's/chr/Uni0/g' \
                | sed 's/psymb/Uni1/g' \
                | sed 's/psyma/Uni2/g' \
                > "contig_order.txt" \
                || { echo "making contig order file for ${SUBSET} failed"; exit 1; }

            ## Make a concatenated reference genome
            ## Note that --gap has to match the setting in the script that
            ## created the dgrp file.
            cp $concat_fasta . \
                || { echo "copying fasta concatenator failed"; exit 1; }
            "./$(basename "$concat_fasta")" \
                --file --gap 1000 "$REFERENCE" "contig_order.txt" | \
                sed 's/[MRWSYKVHDB]/N/g' > \
                "reference.fasta" \
                || { echo "concatenating reference failed"; exit 1; }

            samtools faidx "reference.fasta" \
                || { echo "indexing fasta failed"; exit 1; }

            ## Filter the SNPs to include just the strains of interest
            cp "$filter_dgrp" .
            "./$(basename "$filter_dgrp")" \
                --keep "keep.txt" --min-gt 1 \
                --output "snps.txt" "$SNPCALLS" \
                || { echo "filtering SNPs for failed"; exit 1; }

            ;;

        "run")
            cd "$SCRATCHDIR"
            mkdir -p "$POOL" &&  cd "$POOL"
            echo "Using $(pwd) as tmp for ${POOL}"

            ## Change the reference names and coordinates in the bam file
            ## Note that this may print out a bunch of warnings about reads
            ## with 0 coordinate that samtools is treating as unmapped. These
            ## reads are, in fact, unmapped, so there isn't a problem.
            cp $concat_bam . \
                || { echo "copying bam concatenator failed"; exit 1; }
            samtools view -F 4 -h "${BAMDIR}/${POOL}/alignment.cram" \
                | "./$(basename $concat_bam)" --gap 1000 --file "../contig_order.txt" \
                | samtools view -b - \
                | samtools sort -o "input.bam" - \
                || { echo "concatenating ${POOL} (${BAMDIR}/${POOL}/alignment.cram) failed"; exit 1; }

            samtools index "input.bam" \
                || { echo "indexing ${POOL} bam failed"; exit 1; }

            ## Run HARP
            harp like --bam "input.bam" \
                --region "genome:1-$(cut -f 2 ../reference.fasta.fai)" \
                --refseq "../reference.fasta" \
                --snps "../snps.txt" \
                --stem "genome" \
                || { echo "running likelihood on ${POOL} failed"; exit 1; }

            harp freq --hlk "genome.hlk" --compute_standard_errors \
                --region "genome:1-$(cut -f 2 ../reference.fasta.fai)" \
                || { echo "running freq on ${POOL} failed"; exit 1; }

            rm "genome.hlk" "input.bam" \
                || { echo "removing intermediate files failed"; exit 1; }

            ls . | grep "^samtools" && \
                { echo "found samtools temporary file in ${POOL} output -- failed"; exit 1; }

            cd "$OUTDIR"
            mkdir -p "$POOL" && cd "$POOL"
            rsync -av "${SCRATCHDIR}/${POOL}/" "./" \
                || { echo "copying files from scratch for ${POOL} failed"; exit 1; }

            ;;

        "combine")

            rsync -av "$SCRATCHDIR/" "./" \
                || { echo "copying files from scratch failed"; exit 1; }
            
            ## This will complain about not finding genome.freqs -- this is because
            ## find returns a ".", and ./genome.freqs does not exist. Don't worry
            ## about this error message. Same issue below for genome.freqs.se.
            POOLS=$(find -type 'd' -maxdepth 1 | sed 's/\.\///g' | grep -v "^\.")
            echo "pool" \
                $(head -n 1 "snps.txt" | \
                cut -f 3- -d "," | sed 's/,/ /g') | \
                sed 's/ /\t/g' > "freq.tsv" \
                || { echo "getting header for freqs failed"; exit 1; }
            for p in $POOLS; do
                echo $p $(cut -f 4- -d " " "${p}/genome.freqs") | \
                    sed 's/ /\t/g' >> "freq.tsv" \
                    || { echo "getting freqs from ${p} failed"; exit 1; }
            done
            
            echo "pool" \
                $(head -n 1 "snps.txt" | \
                cut -f 3- -d "," | sed 's/,/ /g') | \
                sed 's/ /\t/g' > "se.tsv" \
                || { echo "getting header for se failed"; exit 1; }
            for p in $POOLS; do
                echo $p $(cut -f 4- -d " " "${p}/genome.freqs.se") | \
                    sed 's/ /\t/g' >> "se.tsv" \
                    || { echo "getting se from ${p} failed"; exit 1; }
            done

            for p in $POOLS; do
                rm -rf "${p}/genome.output/like" \
                    || { echo "deleting like file for ${p} failed"; exit 1; }
            done

            ;;

    esac

    rm "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYDIR"
    mkdir -p "$SCRATCHDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-$(basename $0)-${1}"
    cp "$0" "$sfile"

    i=0
    case $1 in
        "prep")
            echo "$1"$'\t'_ > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        "run")
            for pool in $(ls "$BAMDIR" | grep -v "script\|array\|log\|working"); do
                echo "$1"$'\t'"$pool" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;
        "combine")
            WALLTIME="00:20:00"
            MEM="2GB"
            NTHREADS=1
            echo "$1"$'\t'_ > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        *)
            echo "Unrecognized command"
            ;;
    esac

    if [[ "$i" == 1 ]]; then
        array="0"
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    sbatch --partition="$QUEUE" --time="$WALLTIME" --export=PBS_ENVIRONMENT=PBS_BATCH \
        --mem=$MEM --nodes=1 --ntasks=$NTHREADS --array="$array" \
        "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
