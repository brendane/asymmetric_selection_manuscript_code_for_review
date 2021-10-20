#!/bin/bash
#
# Use BWA mem to align Mtr reads to the Mt5.0 reference genome.
#
# Downloads reads, performs quality trimming and contaminant
# removal before running the alignment.
#
# SETTINGS
#
QUEUE="mangi"
#
# bwa settings - all defaults except MAPQ; placed here for record-keeping
#
THREADS=8           # -t option
MINSEEDLEN=19       # -k option
BANDWIDTH=100       # -w option
ZDROP=100           # -d option
RESEED=1.5          # -r option
DISCARDCOUNT=10000  # -c option
MATCHSCORE=1        # -A option
MMPENALTY=4         # -B option
GAPOPEN=6           # -O option
GAPXPEN=1           # -E option
CLIPPEN=5           # -L option
UNPPEN=9            # -U option
MAPQ=0              # -T option (min alignment score), set to zero to include all alignments

# TrimGalore! settings
# Same settings for all types of reads, regardless of original length
QUALTHRESH=30
LENTHRESH=45
ADAPTSTRINGENCY=3
ERRORRATE=0.1

# bbduk
BBDUK_KMER=25       # Size of kmer to use for contamination filtering
BBDUK_HDIST=0       # Distance between contaminant kmers and read kmers; set to 0 for efficiency
BBDUK_RSKIP=4       # Store only 1/4 of the reference kmers (saves memory)
BBDUK_MEM=24GB


PROJDIR="${HOME}/project/mtr_variant_calls"
REFERENCE="${PROJDIR}/data/assemblies/Mtr_genomes/MtrunA17r5.0-ANR/genome.fasta"
DOWNLOAD_FILE="${PROJDIR}/data/reads/ncbi/sra_runs_to_download_2020-05-20.tsv"
TIM_READS_DIR="${PROJDIR}/data/reads/tim"
CONTAMINANT_DIR="${PROJDIR}/data/assemblies/contaminants"

OUTDIR="${PROJDIR}/results/align/mt5/bwa/2020-05-08"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
ARRAYDIR="${OUTDIR}/arrayjobdata"
SCRIPTDIR="${OUTDIR}/script_copies"
SCRATCHDIR="/scratch.global/${USER}/Mtr_alignment"

check_ncbi_fastq="${PROJDIR}/script/bin/check_ncbi_fastq.py"
de_interleave="${PROJDIR}/script/bin/de_interleave.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bbtools/38.86
    module load bwa/0.7.17
    module load gcc
    module load pigz
    module load samtools/1.9
    module load seqkit/0.13.0
    module load seqtk/1.3
    module load sratoolkit/2.10.7
    module load sratoolkit/2.8.2
    module load trimgalore/0.5.0

    set -euo pipefail

    TASK=$(cut -f 1 "${ARRAYDIR}/${PBS_ARRAYID}")
    ACCESSION=$(cut -f 2 "${ARRAYDIR}/${PBS_ARRAYID}")
    RUN=$(cut -f 3 "${ARRAYDIR}/${PBS_ARRAYID}")
    R1=$(cut -f 4 "${ARRAYDIR}/${PBS_ARRAYID}")
    R2=$(cut -f 5 "${ARRAYDIR}/${PBS_ARRAYID}")

    SECONDS=0

    cd "$SCRATCHDIR"

    ## This won't work if there are any whitespace characters in the
    ## accession numbers, but I know there aren't any.
    runs=($(cut -f 2 "$DOWNLOAD_FILE"))

    case "$TASK" in

        download)
            # Use the SRA toolkit to download reads from NCBI.
            for run in "${runs[@]}"; do
                echo "Downloading ${run}"
                if [[ -f "${run}.fastq.gz" && ! -f "/scratch.global/epste051/ncbi-tmp/sra/${run}.sra" ]]; then
                    echo "    already downloaded, skipping"
                    continue
                fi
                ## Don't know why, but only 2.8.2 was working for me; need a
                ## later version for fasterq-dump
                prefetch.2.8.2 --max-size 100G "$run" \
                    || { echo "fetching ${run} failed"; exit 1; }
                rm -f "${run}.fastq"
                fasterq-dump -e $THREADS --split-spot -O . $run \
                    || { echo "dumping ${run} failed"; exit 1; }
                pigz -p $THREADS --best -f "${run}.fastq" \
                    || { echo "gzipping ${run} failed"; exit 1; }
                rm "/scratch.global/epste051/ncbi-tmp/sra/${run}.sra" \
                    || { echo "deleting cache for ${run} failed"; exit 1; }
            done

            ;;

        split)
            # Change the names and de-interleave; keep runs separate for now.
            cp "$de_interleave" .
            cp "$check_ncbi_fastq" .
            rm -f counts.txt
            for run in "${runs[@]}"; do
                acc="$(awk '$2 == "'"$run"'" { print $1 };' "$DOWNLOAD_FILE")" \
                    || { echo "getting accession name for ${run} failed"; exit 1; }
                echo "$acc, $run"
                seqtk seq -1 "${run}.fastq.gz" > "${acc}-${run}_R1.fastq" \
                    || { echo "extracting 1 reads from ${run} in ${acc} failed"; exit 1; }
                seqtk seq -2 "${run}.fastq.gz" > "${acc}-${run}_R2.fastq" \
                    || { echo "extracting 2 reads from ${run} in ${acc} failed"; exit 1; }
                echo "Finished de-interleaving ${run}"

                ## Count reads pairs
                blwc "${acc}-${run}_R1.fastq" | \
                    sed 's/-/\t/g' | \
                    sed 's/_R1.fastq//g' >> counts.txt \
                    || { echo "counting reads for ${acc} failed"; exit 1; }

                # Check for obvious problems
                "./$("$check_ncbi_fastq")" "${acc}-${run}_R1.fastq" "${acc}-${run}_R2.fastq" \
                    || echo "$acc failed check"

                # gzip
                pigz -p $THREADS -f "${acc}-${run}_R1.fastq" \
                    || { echo "gzipping ${acc}-${run} failed"; exit 1; }
                pigz -p $THREADS -f "${acc}-${run}_R2.fastq" \
                    || { echo "gzipping ${acc}-${run} failed"; exit 1; }
            done

            # Remove temporary files
            mkdir -p srr
            mv SRR* srr/
            ;;

        trim)
            mkdir -p trimmed && cd trimmed
            mkdir -p "$ACCESSION" && cd "$ACCESSION"

            mem="$(Rscript -e 'cat(floor(as.numeric(gsub("[A-Z]", "", "'$BBDUK_MEM'"))*0.85))')g"

            bbduk.sh -Xmx$mem in1="$R1" in2="$R2" \
                out1="${RUN}.1.fastq" out2="${RUN}.2.fastq" \
                stats="bbduk_stats.${RUN}.txt" \
                ref="${CONTAMINANT_DIR}/human/genome.fasta,${CONTAMINANT_DIR}/ensifer/genome.fasta,${CONTAMINANT_DIR}/phix/genome.fasta" \
                k=$BBDUK_KMER hdist=$BBDUK_HDIST rskip=$BBDUK_RSKIP \
                || { echo "bbduk on ${ACCESSION}, ${RUN} failed"; exit 1; }

            o1="${RUN}.1_val_1.fq"
            o2="${RUN}.2_val_2.fq"
            trim_galore --quality "$QUALTHRESH" \
                --stringency "$ADAPTSTRINGENCY" \
                -e "$ERRORRATE" \
                --length "$LENTHRESH" --paired \
                "${RUN}.1.fastq" "${RUN}.2.fastq" | tee "log.${RUN}.txt" \
                || { echo "running trim galore on ${ACCESSION}, ${RUN} failed"; exit 1; }

            mv "$o1" "paired.${RUN}.r1.fastq" \
                || { echo "moving R1 output for ${ACCESSION}, ${RUN} failed"; exit 1; }
            mv "$o2" "paired.${RUN}.r2.fastq" \
                || { echo "moving R1 output for ${ACCESSION}, ${RUN} failed"; exit 1; }
            rm "${RUN}.1.fastq" "${RUN}.2.fastq"
            ;;

        align)
            cd trimmed
            cd "$ACCESSION"
            mkdir -p "${OUTDIR}/${ACCESSION}"

            # Align the paired end reads.
            r1="paired.${RUN}.r1.fastq"
            r2="paired.${RUN}.r2.fastq"
            rg="@RG\\tID:${RUN}\\tSM:${ACCESSION}"
            if [[ "$ACCESSION" =~ "A_PI" ]]; then
                rg="@RG\\tID:${ACCESSION}.${RUN}\\tSM:${ACCESSION}"
            fi
            bwa mem -R "$rg" \
                -t "$THREADS" -k "$MINSEEDLEN" -w "$BANDWIDTH" -d "$ZDROP" \
                -r "$RESEED" -c "$DISCARDCOUNT" -A "$MATCHSCORE" -B "$MMPENALTY" \
                -O "$GAPOPEN" -E "$GAPXPEN" -L "$CLIPPEN" -U "$UNPPEN" -T "$MAPQ" \
                "$REFERENCE" "$r1" "$r2" \
                | samtools view -b \
                | samtools sort -n \
                > "initial.${RUN}.bam" \
                || { echo "aligning ${RUN} ${ACCESSION} failed"; exit 1; }

            ## Mark duplicates, sort by position, and convert to CRAM
            samtools fixmate --reference "$REFERENCE" -m "initial.${RUN}.bam" - \
                | samtools sort - \
                | samtools markdup --reference "$REFERENCE" - - \
                | samtools sort -O "CRAM" --reference "$REFERENCE" - \
                > "alignment.${RUN}.cram" \
                || { echo "marking duplications in ${RUN}, ${ACCESSION} failed"; exit 1; }

            # Index the cram file
            samtools index "alignment.${RUN}.cram" \
                || { echo "indexing ${RUN} ${ACCESSION} failed"; exit 1; }

            # Various coverage and read length statistics
            samtools stats "alignment.${RUN}.cram" > "stats.${RUN}.txt" \
                || { echo "stats failed on ${RUN} ${ACCESSION}"; exit 1; }

            # Check for truncated file or invalid format -- just in case of
            # uncaught errors
            samtools flagstat "alignment.${RUN}.cram" | \
                grep -i "invalid\\|truncat" &&
                { echo "invalid cram file for ${RUN} ${ACCESSION} failed"; exit 1; }

            cp "alignment.${RUN}.cram" "alignment.${RUN}.cram.crai" \
                "stats.${RUN}.txt" "${OUTDIR}/${ACCESSION}/" \
                || { echo "copying files from scratch for ${RUN} ${ACCESSION} failed"; exit 1; }

            ;;

        stats)
            seqkit stats -T -j 24 *.fastq.gz "$TIM_READS_DIR/"*.fastq.gz \
                > "${OUTDIR}/raw_read_counts.tsv" \
                || { echo "counting raw reads failed"; exit 1; }

            cd trimmed
            grep "Matched" */bbduk_stats.*.txt \
                | sed 's/\/bbduk_stats\./\t/g' \
                | sed 's/.txt:#Matched//g' \
                > "${OUTDIR}/bbduk_contaminant_stats.txt" \
                || { echo "making table of bbduk contaminant percentage failed"; exit 1; }
            cd ..

            seqkit stats -T -j 24 trimmed/*/*.fastq \
                > "${OUTDIR}/trimmed_read_counts.tsv" \
                || { echo "counting trimmed reads failed"; exit 1; }

            cd trimmed
            grep "reads unmapped" */stats.*.txt \
                | sed 's/\/stats./\t/g' \
                | sed 's/.txt:SN\treads unmapped:\t/\t/g' \
                > "${OUTDIR}/unmapped_read_counts.tsv" \
                || { echo "making table of unmapped reads failed"; exit 1; }
            grep "reads mapped:" */stats.*.txt \
                | sed 's/\/stats./\t/g' \
                | sed 's/.txt:SN\treads unmapped:\t/\t/g' \
                > "${OUTDIR}/mapped_read_counts.tsv" \
                || { echo "making table of mapped reads failed"; exit 1; }
            cd ..
            ;;

    esac

    echo "RUN TIME $SECONDS ($(($SECONDS/60)) minutes) for ${ACCESSION}"

    rm "${ARRAYDIR}/${PBS_ARRAYID}"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYDIR"
    mkdir -p "$SCRATCHDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-${1}-$(basename $0)"
    cp "$0" "$sfile"

    i=0

    MEM=4GB

    case "$1" in
        download)
            WALLTIME=48:00:00
            echo "$1" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;

        split)
            WALLTIME=72:00:00
            echo "$1" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;

        trim)
            WALLTIME=10:00:00
            MEM=$BBDUK_MEM
            for r1file in "${SCRATCHDIR}/"*_R1.fastq.gz; do 
                acc="$(basename "$r1file" | sed 's/-.\+//g')"
                run="$(basename "$r1file" | sed 's/.\+-//g' | sed 's/_.\+//g')"
                r2file="${r1file/_R1/_R2}"
                echo "$1"$'\t'"$acc"$'\t'"$run"$'\t'"$r1file"$'\t'"$r2file" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            for r1file in "${TIM_READS_DIR}"/*_R1.fastq.gz; do
                acc="$(basename "$r1file" | sed 's/20170317.//g' | sed 's/_R1.fastq.gz//g' | sed 's/-/_/g')"
                run=20170317
                r2file="${r1file/_R1/_R2}"
                echo "$1"$'\t'"$acc"$'\t'"$run"$'\t'"$r1file"$'\t'"$r2file" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        align)
            WALLTIME=24:00:00
            MEM=48GB
            for accdir in "${SCRATCHDIR}/trimmed"/*; do
                acc="$(basename "$accdir")"
                for r1file in "$accdir"/*.r1.fastq; do
                    run="$(basename "${r1file/.r1.fastq/}" | sed 's/^paired.//g')"
                    echo "$1"$'\t'"$acc"$'\t'"$run"$'\t'"_"$'\t'"_" > "${ARRAYDIR}/${i}"
                    i=$(($i+1))
                done
            done
            ;;

        stats)
            MEM=10GB
            THREADS=24
            WALLTIME=02:30:00
            echo "$1" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
    esac

    if [[ "$i" == 1 ]]; then
        array="0"
    else
        array="0-$(($i-1))"
    fi
    array="$(ls "$ARRAYDIR" | tr '\n' ',' | sed 's/,$//g')"

    cd "$WORKDIR"
    qsub2 -A "tiffinp" -W group_list="tiffinp" \
        -q "$QUEUE" -l "walltime=$WALLTIME,nodes=1:ppn=$THREADS,mem=$MEM" -t "$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
