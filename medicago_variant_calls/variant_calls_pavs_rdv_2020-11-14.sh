#!/bin/bash
#
# Use a simple cutoff based read depth variation method to detect CNVs
# relative to the Mt5.0 reference genome. Based on the method
# described in Joseph Guhlin's thesis.
#
# Uses only libraries sequenced on Illumina GA IIx machines because
# these early libraries had different coverage characteristics than
# later libraries, and they are the most common library in the dataset.
#
# Normalization: There is no need to divide by gene length because each
# gene is only compared to itself. FPKM usually uses the total number of
# mapped reads as the denominator, so I counted up the reads that mapped to
# a gene model, and used that for the library size.
#
# SETTINGS
#
QUEUE=mangi
#
MIN_MAPQ=30
MODE=union
FEATURE_ID=locus_tag
FEATURE_TYPE=gene
STRANDED=no
NONUNIQUE=none
#

tsv2vcf="""
import csv
import sys
annot = {}
with open(sys.argv[1], 'rt') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for row in rdr:
        if row[0].startswith('#'): continue
        if row[2] != 'gene': continue
        chrom = row[0]; pos = row[3]
        lt = {x.split('=')[0]:x.split('=')[1] for x in row[8].split(';')}['locus_tag']
        annot[lt] = (chrom, pos)
with open(sys.argv[2], 'rt') as handle:
    genotypes = handle.readline().strip().split('\t')[1:]
    print('##fileformat=VCFv4.3')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">')
    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(genotypes))
    rdr = csv.reader(handle, delimiter='\t')
    for row in rdr:
        chrom, pos = annot[row[0]]
        genotypes = []
        up = 0
        down = 0
        norm = 0
        for r in row[1:]:
            try:
                rr = float(r)
            except:
                genotypes.append('N')
                continue
            if -1 >= rr > -2 or 1 <= rr < 2:
                genotypes.append('N')
            elif rr <= -2:
                genotypes.append('A')
                down += 1
            elif rr >= 2:
                genotypes.append('G')
                up += 1
            else:
                genotypes.append('T')
                norm += 1
        output_genotypes = []
        geno_code = {'T':'0', 'N':'.'}
        alts = []
        if up > 0:
            geno_code['G'] = '1'; geno_code['A'] = '2'
            alts.append('G')
            if down > 0: alts.append('A')
        else:
            geno_code['A'] = '1'
            alts.append('A')
        for g in genotypes:
            output_genotypes.append(geno_code[g])
        print(chrom +'\t' + pos + '\t' + row[0] + '\tT\t' + ','.join(alts)
              + '\t.\t.\t.\tGT\t' + '\t'.join(output_genotypes))
"""

PROJDIR="${HOME}/project/mtr_variant_calls"
TIER2_ALN_DIR="s3://medicago_variant_calls/trim_and_align/2020-05-08"
ACCESSIONS_TABLE="${PROJDIR}/data/reads/ncbi/medicago_reads_ncbi.csv"
OUTDIR="${PROJDIR}/results/variant_calls/pavs/rdv/2020-11-14"
REFDIR="${PROJDIR}/data/assemblies/Mtr_genomes/MtrunA17r5.0-ANR"
REFERENCE="${REFDIR}/genome.fasta"
ANNOTATIONS="${REFDIR}/genome.gff"

WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYJOBDIR="${OUTDIR}/arraydata"
SCRATCHDIR="/scratch.global/${USER}/variant_calls_pavs_MTR_rdv_2020-11-14"

if [[ "$PBS_ENVIRONMENT" == PBS_BATCH ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bcftools/1.9
    module load htseq/0.11.0
    module load R/3.6.3
    module load samtools/1.10

    SECONDS=0

    TASK="$(cut -f 1 "${ARRAYJOBDIR}/${PBS_ARRAYID}")"
    GENOTYPE="$(cut -f 2 "${ARRAYJOBDIR}/${PBS_ARRAYID}")"

    cd "$SCRATCHDIR"

    case "$TASK" in

        "prep")
            ## List the genotypes currently in tier 2 storage
            s3cmd ls "$TIER2_ALN_DIR/" \
                | grep "HM\|A_PI" \
                | sed 's/\/$//g' \
                | sed 's/.\+\///g' \
                > genotypes_with_aln.txt \
                || { echo "getting list of available genotypes failed"; exit 1; }

            ## Pick out only Illumina GA IIx runs for consistency across
            ## genotypes
            grep "IIx" "$ACCESSIONS_TABLE" \
                | cut -f 1 -d, \
                > gaIIx_accessions.txt \
                || { echo "getting list of SRA accessions with matching technology failed"; exit 1; }

            ## Download reads from tier 2
            while read -r genotype; do
                mkdir -p "$genotype"
                echo "$genotype"
                s3cmd ls "${TIER2_ALN_DIR}/${genotype}/" \
                    | grep "cram$"  \
                    | sed 's/.\+\.SRR/SRR/g' \
                    | sed 's/\.cram$//g' \
                    > srr \
                    || { echo "no alignments for ${genotype}, skipping"; rm -r "$genotype"; continue; }
                while read -r sra; do
                    echo "    $sra"
                    grep "$sra" gaIIx_accessions.txt \
                        && s3cmd get --force "${TIER2_ALN_DIR}/${genotype}/alignment.${sra}.cram" \
                        "$genotype/" \
                        || { echo "${sra} not among Illumina GA IIx runs, skipping"; continue; }
                done < srr
            done < genotypes_with_aln.txt
            rm srr
            ;;

        "coverage")
            cd "$GENOTYPE"

            ## Merge bams
            samtools merge -u -O BAM --reference "$REFERENCE" \
                - alignment.*.cram \
                | samtools sort --reference "$REFERENCE" -n -l 0 - \
                > merged.bam \
                || { echo "merging alignments for ${GENOTYPE} failed"; exit 1; }

            ## Calculate coverage for each gene
            awk -F$'\t' '$3 == "gene" { print $0 };' "$ANNOTATIONS" \
                > genes.gff \
                || { echo "extracting genes failed"; exit 1; }
            PYTHONPATH="" htseq-count -f bam -r name \
                --nonunique $NONUNIQUE \
                -s $STRANDED \
                -a $MIN_MAPQ \
                -t $FEATURE_TYPE \
                -i $FEATURE_ID \
                -m $MODE \
                merged.bam genes.gff \
                > "count.txt" \
               || { echo "counting coverage for ${GENOTYPE} failed"; exit 1; }

            ## Check for failed samtools runs (note could be from previous
            ## runs of the script)
            grep "^samtools" && { echo "found samtools temp files in ${GENOTYPE}, failed"; exit 1; }

            rm genes.gff
            ;;

        "combine")
#            ## Calculate coverage relative to HM101 for every gene and genotype
#            ## Need to divide by total number of reads but not by gene
#            ## length b/c we are not comparing genes to each other.
#            for genotype in HM*; do
#                total=$(grep -v "^__" "${genotype}/count.txt" | awk -F$'\t' 'BEGIN {x=0}; {x+=$2}; END {print x};')
#                grep -v "^__" "${genotype}/count.txt" \
#                    | awk -F$'\t' 'BEGIN {print "gene\t'$genotype'"}; {print $1 "\t" $2 * 1000000 / '$total'};' \
#                    > "${genotype}/normalized.tsv" \
#                    || { echo "normalizing ${genotype} by total read count failed"; exit 1; }
#            done
#
#            paste */normalized.tsv > normalized_by_count.tsv \
#                || { echo "putting all count normalized files together failed"; exit 1; }

            div_by_hm101=''
            div_by_hm101+='x = read.csv("normalized_by_count.tsv", sep="\t", header=TRUE, as.is=TRUE); '
            div_by_hm101+='x = x[, !(colSums(is.na(x))==nrow(x))]; '
            div_by_hm101+='for(i in seq(1, ncol(x), 2)) {stopifnot(all(x[, 1] == x[, i]))}; '
            div_by_hm101+='y = x[, c(1, seq(2, ncol(x), 2))]; '
            div_by_hm101+='z = apply(y[, -1], 2, function(r) { rr=log2(r / y[, "HM101"]); '
            div_by_hm101+='ifelse(y[, "HM101"]==0 & is.infinite(rr), 10, ifelse(is.infinite(rr) & !(is.na(rr)) & r==0, -10, r))}); '
            div_by_hm101+='write.table(cbind(gene=y[, "gene"], z), stdout(), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE); '
            Rscript -e "$div_by_hm101" > normalized_by_HM101.tsv \
                || { echo "normalizing by HM101 failed"; exit 1; }
            ;;


        vcf)
            ## Convert to genotype calls and VCF format
            python -c "$tsv2vcf" "$ANNOTATIONS" normalized_by_HM101.tsv \
                > pavs.vcf \
                || { echo "calling PAVs failed"; exit 1; }

            cp normalized_by* pavs.vcf "$OUTDIR"
            ;;

    esac

    rm "${ARRAYJOBDIR}/${PBS_ARRAYID}"

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

    i=0
    case "$1" in
        "prep")
            NTHREADS=1
            WALLTIME=06:00:00
            MEM=2GB
            echo "$1"$'\t'"_" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;
        "coverage")
            WALLTIME=03:00:00
            MEM=8GB
            NTHREADS=1
            for genotype in $(ls "$SCRATCHDIR" | grep "HM"); do
                echo "$1"$'\t'"$genotype" > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            ;;
        "combine")
            NTHREADS=1
            WALLTIME=00:10:00
            MEM=4GB
            echo "$1"$'\t'"_" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;
        "vcf")
            NTHREADS=1
            WALLTIME=00:05:00
            MEM=4GB
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
    qsub2 -A "tiffinp" -W group_list="tiffinp" -q "$QUEUE" \
        -l "walltime=${WALLTIME},nodes=1:ppn=${NTHREADS},mem=${MEM}" \
        -t "$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"
fi
