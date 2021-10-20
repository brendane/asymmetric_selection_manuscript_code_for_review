#!/bin/bash
#
# Project allele frequencies using relative fitness, assuming that the
# strains were at the same frequency at the beginning of the experiment.
#
# UPDATE 2021-05-27: Removed allele frequency code and re-ran with
# updated strain frequencies
#
# SETTINGS
#
# pbs
QUEUE="amdsmall"
#

project_af=''
project_af+='argv = commandArgs(trailingOnly=TRUE); '
project_af+='infile = argv[1]; '
project_af+='community = basename(gsub(".tsv$", "", infile)); '
project_af+='X = read.csv(infile, sep="\t", header=TRUE, as.is=TRUE, check.names=FALSE); '
project_af+='envs = unique(gsub("_Fit_.+", "", colnames(X)[grepl("_Fit", colnames(X))])); '
project_af+='if(length(envs) == 0) { '
project_af+='cols = colnames(X)[!grepl("minus|^strain|_freq$|_var$|_fitdiff$|_fitbin$", colnames(X))]; '
project_af+='} else { '
project_af+='cols = paste0(envs, "_Fit_med"); '
project_af+='use_mean = !(cols %in% colnames(X)); '
project_af+='cols[use_mean] = paste0(envs[use_mean], "_Fit_mean"); '
project_af+='}; '
project_af+='med_fit = as.matrix(X[, cols]); '
project_af+='colnames(med_fit) = gsub("_Fit_.+", "", colnames(med_fit)); '
project_af+='colnames(med_fit) = gsub("_fit$", "", colnames(med_fit)); '
project_af+='rownames(med_fit) = X[, "strain"]; '
project_af+='proj_freq = apply(med_fit, 2, function(x) { (1/length(x)) * (2^x) }); '
project_af+='proj_freq = proj_freq[, apply(proj_freq, 2, function(x) sum(!is.na(x))) == nrow(proj_freq)]; '
project_af+='write.table(cbind("strain"=rownames(proj_freq), proj_freq), file=stdout(), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE); '

PROJDIR="${HOME}/project/select_reseq"
REL_FIT_DIR="${PROJDIR}/results/relative_fitness/compile/2020-05-23"

OUTDIR="${PROJDIR}/results/phenotypes/rhizobia_phenotypes/projected_frequency/2021-01-02"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYJOBDIR="${OUTDIR}/arrayjobdata"
LOGDIR="${OUTDIR}/log"
WORKDIR="${OUTDIR}/working"
SCRATCHDIR="/scratch.global/${USER}/projected_af_2021-01-02"

ld="${PROJDIR}/script/bin/intergenomic_ld"
eucl="${PROJDIR}/script/bin/euclidean_genetic_distance_af.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bcftools/1.10.2
    module load htslib/1.9
    module load R/3.4.3

    set -euo pipefail

    datafile="${ARRAYJOBDIR}/${SLURM_ARRAY_TASK_ID}"
    TASK="$(cut -f 1 "$datafile")"
    SPECIES="$(cut -f 2 "$datafile")"
    COMMUNITY="$(cut -f 3 "$datafile")"

    SECONDS=0

    cd "$SCRATCHDIR"
    mkdir -p "$SPECIES" && cd "$SPECIES"
    mkdir -p "$OUTDIR/$SPECIES"

    case "$TASK" in

        allele-frequency)

            ## If community == SP1, extract just the strains that belong
            ## to this species
            if [[ "$COMMUNITY" == "SP1" ]]; then
                grep "^strain\|${SPECIES}" "${REL_FIT_DIR}/${COMMUNITY}.tsv" \
                    > relfit \
                    || { echo "filtering relative fitness for ${COMMUNITY}, ${SPECIES} failed"; exit 1; }
            else
                cp "${REL_FIT_DIR}/${COMMUNITY}.tsv" relfit \
                    || { echo "copying relative fitness for ${COMMUNITY}, ${SPECIES} failed"; exit 1; }
            fi

            # Calculate strain frequences
            Rscript -e "$project_af" relfit \
                > "${COMMUNITY}.projected_strain_frequencies.tsv" \
                || { echo "calculating rhizobial strain frequencies for ${COMMUNITY} failed"; exit 1; }

            cp "$COMMUNITY".* "${OUTDIR}/${SPECIES}"

    esac

    echo "RUN TIME = ${SECONDS} seconds or $(($SECONDS/60)) minutes"

    rm -f "$datafile"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYJOBDIR"
    mkdir -p "$SCRATCHDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-${1}-$(basename $0)"
    cp $0 $sfile

    i=0
    case "$1" in

        allele-frequency)
            WALLTIME=00:05:00
            NTHREADS=1
            MEM=1G
            for species in medicae meliloti; do
                if [[ "$species" == "medicae" ]]; then
                    commfile="$REL_FIT_DIR/SP1.tsv"
                    comm="$(basename "${commfile/.tsv/}")"
                    echo "$1"$'\t'"$species"$'\t'"$comm" > "${ARRAYJOBDIR}/${i}"
                    i=$(($i+1))
                else
                    for commfile in "$REL_FIT_DIR/"*.tsv; do
                        comm="$(basename "${commfile/.tsv/}")"
                        echo "$1"$'\t'"$species"$'\t'"$comm" > "${ARRAYJOBDIR}/${i}"
                        i=$(($i+1))
                    done
                fi
            done
            ;;

    esac

    if [[ "$i" == 1 ]]; then
        array=0
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    sbatch --account=tiffinp --export=PBS_ENVIRONMENT=PBS_BATCH \
        --partition=$QUEUE --nodes=1 --ntasks=$NTHREADS \
        --time=$WALLTIME --mem=$MEM --array="$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
