#!/bin/bash
#
# Trim reads from the SpanFran2 winter 2020 GWAS experiment.
#
# SETTINGS
#
QUEUE="small"
WALLTIME="03:30:00"
#
QUALTYPE="sanger"
QUALTHRESH=30
LENTHRESH=99
SLENTHRESH=100
ADAPTSTRINGENCY=3
ERRORRATE=0.1

PROJDIR="${HOME}/project/select_reseq"
READSDIR="${PROJDIR}/data/reads/SpanFran2/gwas"
OUTDIR="${PROJDIR}/results/trim_reads/spanfran2/trim_galore/2020-10-23_gwas"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
ARRAYDIR="${OUTDIR}/arrayjobdata"
SCRIPTDIR="${OUTDIR}/script_copies"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load gcc
    module load htslib/1.9
    module load trimgalore/0.5.0

    set -euo pipefail

    SAMPLE="$(cut -f 1 -d " " "${ARRAYDIR}/${PBS_ARRAYID}")"

    cd "$OUTDIR"
    mkdir -p "$SAMPLE"
    cd "$SAMPLE"

    R1="${READSDIR}/${SAMPLE}_R1.fastq.gz"
    R2="${READSDIR}/${SAMPLE}_R2.fastq.gz"
    o1=$(basename $R1 | sed 's/\\.fastq\\.gz/_val_1.fq.gz/g')
    o2=$(basename $R2 | sed 's/\\.fastq\\.gz/_val_2.fq.gz/g')
    u1=$(basename $R1 | sed 's/\\.fastq\\.gz/_unpaired_1.fq.gz/g')
    u2=$(basename $R2 | sed 's/\\.fastq\\.gz/_unpaired_2.fq.gz/g')

    # Trim with TrimGalore!
    trim_galore --quality $QUALTHRESH --stringency $ADAPTSTRINGENCY \
        -e $ERRORRATE --gzip --length $LENTHRESH --paired \
        --retain_unpaired --length_1 $SLENTHRESH --length_2 $SLENTHRESH \
        $R1 $R2 | tee "log.txt" \
        || { echo "running trim galore on ${SAMPLE} failed"; exit 1; }

    echo "cutadapt version:" $(cutadapt --version) >> "log.txt"

    zcat "$u1" "$u2" | bgzip -c > \
        "single.fastq.gz" \
        || { echo "concatenating singles for ${SAMPLE} failed"; exit 1; }
    rm "$u1" "$u2"

    mv "$o1" "paired.r1.fastq.gz" \
        || { echo "moving R1 output for ${SAMPLE} failed"; exit 1; }
    mv "$o2" "paired.r2.fastq.gz" \
        || { echo "moving R1 output for ${SAMPLE} failed"; exit 1; }

    # Count reads
    rm -f "read_counts.txt"
    echo "pairs	$(zcat "paired.r1.fastq.gz" | blwc - | sed 's/pairs.\\+\\t/pairs/g')" > "read_counts.txt" \
        || { echo "counting number of read pairs left in ${STRAIN} failed"; exit 1; }
    echo "singles	$(zcat "single.fastq.gz" | blwc - | sed 's/singles.\\+\\t/singles/g')" >> "read_counts.txt" \
        || { echo "counting number of single read left in ${STRAIN} failed"; exit 1; }

    rm "${ARRAYDIR}/${PBS_ARRAYID}"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-$(basename $0)"
    cp "$0" "$sfile"

    i=0
    ## Note for future: if there are filenames with whitespace, this
    ## would not work properly.
    for sample in $(ls "$READSDIR" | grep "_R1.fastq.gz" | sed 's/_R1.fastq.gz//g'); do
        echo "$sample" > "${ARRAYDIR}/${i}"
        i=$(($i+1))
    done

    if [[ "$i" == 1 ]]; then
        array="0"
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    qsub2 -A "tiffinp" -W group_list="tiffinp" \
        -q "$QUEUE" -l walltime="$WALLTIME" -t "$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
