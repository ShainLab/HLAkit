#!/usr/bin/env bash

set -euo pipefail

# run mpileup

usage() {
    cat << EOF
Usage: $0 --resultdir <dir> --resultdir <path> --somatic_mutations <file> --fastafile <path>

Required arguments:
    -o, --resultdir <dir>           Directory with RNA-Seq bam files
    -m, --somatic_mutations <file>  Somatic mutations file
    -f, --fastafile <path>          Path to fasta file. fasta.fai should be present in the same directory
    -h, --hlakit <path>             Path to hlakit

Optional arguements:
    -h, --help                      Show this help message


Example:
    $0 --resultdir /path/to/resultdirectory --somatic_mutations somatic_mutations.txt --fastafile /path/to/fastafile
    $0 -o /path/to/resultdirectory -m somatic_mutations.txt -f /path/to/fastafile
EOF
    exit 1
}

# Initialize variables
resultdir=""
somatic_mutations=""
fastafile=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--resultdir)
            resultdir="$2"
            shift 2
            ;;
        -m|--somatic_mutations)
            somatic_mutations="$2"
            shift 2
            ;;
        -f|--fastafile)
            fastafile="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [ -z "$resultdir" ] || [ -z "$somatic_mutations" ] || [ -z $fastafile ]; then
    echo "Error: All arguments are required"
    usage
fi

# Validate resultdir exists
if [ ! -d "$resultdir" ]; then
    echo "Error: Result directory does not exist: $resultdir"
    exit 1
fi

# Validate somatic mutations file
if [ ! -f "$somatic_mutations" ]; then
    echo "Error: Somatic mutations file not found: $somatic_mutations"
    exit 1
elif [ ! -s "$somatic_mutations" ]; then
    echo "Error: Somatic mutations file is empty: $somatic_mutations"
    exit 1
fi

# Validate fastafile
if [ ! -f "$fastafile" ]; then
    echo "Error: FASTA file not found at expected path: $fastafile"
    exit 1
fi
if [ ! -f "${fastafile}.fai" ]; then
    echo "Error: FASTA index not found: ${fastafile}.fai — run 'samtools faidx $fastafile' before running this script."
    exit 1
fi

RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m'
log_warn() {
    local msg="[$(date +'%Y-%m-%d %H:%M:%S')] [WARNING] $*"
    echo -e "${YELLOW}${msg}${NC}" >&2
}


# Validate samtools is available
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found in PATH. Please load the samtools module or activate the correct environment."
    exit 1
fi

rna_tumor_mapq0_mpileupout="$resultdir/rna_tumor_mapq0_mpileupout.txt"
rna_tumor_mapqnonzero_mpileupout="$resultdir/rna_tumor_mapqnonzero_mpileupout.txt"


:>$rna_tumor_mapq0_mpileupout
:>$rna_tumor_mapqnonzero_mpileupout

#if no newline char at the end of somatic mutations review file, add \n, else while loop will skip the last variant
if [ -n "$(tail -c 1 $somatic_mutations)" ]
then
    echo >> $somatic_mutations
fi

# --- Indel position correction ---------------------------------------------
mpileup_input="$somatic_mutations"

if awk -F'\t' 'NR > 1 && NF >= 5 && (length($4) != length($5) || $4 == "-" || $4 == "." || $5 == "-" || $5 == ".") {f = 1; exit} END {exit !f}' "$somatic_mutations"; then
    indel_bam="$resultdir/Aligned.sortedByCoord.out_mapqcorrected.bam"
    if [ ! -f "$indel_bam" ]; then
        echo "Error: indels present but BAM for position correction not found: $indel_bam"
        exit 1
    fi
    if ! command -v Rscript &> /dev/null; then
        echo "Error: Rscript not found in PATH (needed for indel position correction)."
        exit 1
    fi

    mpileup_input="$resultdir/somatic_mutations.rna_indelpos.txt"
    echo "Indels found — finding RNA-supported indel positions ..."
    Rscript "$hlakit/rna_indel_pos_correction.R" \
        --mutations "$somatic_mutations" \
        --bam "$indel_bam" \
        --outfile "$mpileup_input"

    if [ ! -s "$mpileup_input" ]; then
        echo "Error: indel position correction produced no output: $mpileup_input"
        exit 1
    fi
fi
# -----------------------------------------------------------------------------

echo "Running mpileup to count ref and alt bases of mutations ..."

n_processed=0
n_missing_bam=0

# run mpileup
while read LINE
do
    if echo $LINE | egrep -q -v "CHROM"
    then
        read -r chr pos id ref alt rest <<< "$LINE"
        region=$(echo $chr":"$pos"-"$pos)

        rna_tumor_mapq0_bam=$resultdir/rna.tumor.MAPQzero.bam
        rna_tumor_mapqnonzero_bam=$resultdir/rna.tumor.MAPQnonzero.bam

        reflen=${#ref}
        altlen=${#alt}

        ref=${ref:0:1}
        alt=${alt:0:1}

        # Validate BAM files exist before running mpileup
        bam_ok=1
        for bam in "$rna_tumor_mapq0_bam" "$rna_tumor_mapqnonzero_bam"; do
            if [ ! -f "$bam" ]; then
                log_warn "BAM file not found $bam."
                bam_ok=0
                n_missing_bam=$((n_missing_bam + 1))
                break
            fi
            if [ ! -f "${bam}.bai" ] && [ ! -f "${bam%.bam}.bai" ]; then
                log_warn "BAM index (.bai) not found for $bam — samtools mpileup may fail."
            fi
        done
        if [ "$bam_ok" -eq 0 ]; then
            continue
        fi

        echo "Running mpileup on allele: $chr pos: $pos ..."
        if [ $reflen -lt $altlen ]; then
            samtools mpileup -A -Q 0 -r $region -f $fastafile -aa --output-QNAME $rna_tumor_mapq0_bam | cut -f1-5,7- >> $rna_tumor_mapq0_mpileupout
            samtools mpileup -A -Q 0 -r $region -f $fastafile -aa --output-QNAME $rna_tumor_mapqnonzero_bam | cut -f1-5,7- >> $rna_tumor_mapqnonzero_mpileupout
        else
            samtools mpileup -r $region -f $fastafile -aa --output-QNAME $rna_tumor_mapq0_bam | cut -f1-5,7- >> $rna_tumor_mapq0_mpileupout
            samtools mpileup -r $region -f $fastafile -aa --output-QNAME $rna_tumor_mapqnonzero_bam | cut -f1-5,7- >> $rna_tumor_mapqnonzero_mpileupout
        fi
        n_processed=$((n_processed + 1))
    fi

done < "$mpileup_input"

echo "Mpileup complete: $n_processed variant(s) processed."
if [ "$n_missing_bam" -gt 0 ]; then
    log_warn "$n_missing_bam variant(s) skipped due to missing BAM files."
fi

if [ "$n_processed" -eq 0 ]; then
    echo "Error: No variants were processed. Check that the somatic mutations file has data rows beyond the header, and that BAM files exist in $resultdir."
    exit 1
fi

# Validate mpileup outputs are non-empty before awk post-processing
for f in "$rna_tumor_mapq0_mpileupout" "$rna_tumor_mapqnonzero_mpileupout"; do
    if [ ! -s "$f" ]; then
        echo "Error: Mpileup output file is empty after processing: $f — check BAM files, region coordinates, and FASTA file."
        exit 1
    fi
done

awk '{gsub(/\$/, "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $rna_tumor_mapq0_mpileupout > ${rna_tumor_mapq0_mpileupout/.txt/.tmp.txt}
mv ${rna_tumor_mapq0_mpileupout/.txt/.tmp.txt} $rna_tumor_mapq0_mpileupout
awk '{gsub(/\$/, "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $rna_tumor_mapqnonzero_mpileupout > ${rna_tumor_mapqnonzero_mpileupout/.txt/.tmp.txt}
mv ${rna_tumor_mapqnonzero_mpileupout/.txt/.tmp.txt} $rna_tumor_mapqnonzero_mpileupout


echo "Done!"
