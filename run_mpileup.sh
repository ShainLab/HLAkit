#!/bin/bash

set -euo pipefail

# run mpileup

usage() {
    cat << EOF
Usage: $0 --resultdir <dir> --threads <int>

Required arguments:
    -o, --resultdir <dir>           Output directory for results
    -m, --somatic_mutations <file>  Somatic mutations file
    -p, --hlakit <path>        Path to hlakit

Optional arguements:
    -h, --help                Show this help message


Example:
    $0 --resultdir /path/to/resultdirectory --somatic_mutations somatic_mutations.txt --hlakit /path/to/hlakit
    $0 -o /path/to/resultdirectory -m somatic_mutations.txt -p /path/to/hlakit
EOF
    exit 1
}

# Initialize variables
resultdir=""
somatic_mutations=""
hlakit=""

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
        -p|--hlakit)
            hlakit="$2"
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
if [ -z "$resultdir" ] || [ -z "$somatic_mutations" ] || [ -z $hlakit ] ; then
    echo "Error: All arguments are required"
    usage
fi

RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m'
log_warn() {
    local msg="[$(date +'%Y-%m-%d %H:%M:%S')] [WARNING] $*"
    echo -e "${YELLOW}${msg}${NC}" >&2
}



fastafile="$hlakit/resources/hla.fasta"
normal_mpileupout="$resultdir/normal_mpileupout.txt"
tumor_mapq0_mpileupout="$resultdir/tumor_mapq0_mpileupout.txt"
tumor_mapqnonzero_mpileupout="$resultdir/tumor_mapqnonzero_mpileupout.txt"


:>$normal_mpileupout
:>$tumor_mapq0_mpileupout
:>$tumor_mapqnonzero_mpileupout

#if no newline char at the end of somatic mutations review file, add \n, else while loop will skip the last variant
if [ -n "$(tail -c 1 $somatic_mutations)" ]
then
    echo >> $somatic_mutations
fi

echo "Running mpileup to count ref and alt bases of mutations ..."

# run mpileup
while read LINE
do
    if echo $LINE | egrep -q -v "CHROM"
    then
        read -r chr pos id ref alt rest <<< "$LINE"
        region=$(echo $chr":"$pos"-"$pos)

        normal_snp_bam=$resultdir/${chr}.normal.SNP.bam
        tumor_mapq0_snp_bam=$resultdir/${chr}.tumor.SNP.MAPQzero.bam
        tumor_mapqnonzero_snp_bam=$resultdir/${chr}.tumor.SNP.MAPQnonzero.bam
        normal_dnp_bam=$resultdir/${chr}.normal.DNP.bam
        tumor_mapq0_dnp_bam=$resultdir/${chr}.tumor.DNP.MAPQzero.bam
        tumor_mapqnonzero_dnp_bam=$resultdir/${chr}.tumor.DNP.MAPQnonzero.bam

        reflen=${#ref}
        altlen=${#alt}

        ref=${ref:0:1}
        alt=${alt:0:1}
        
        if [ $reflen -eq $altlen ] && [ $reflen -eq 2 ]; then
            normal_bam=$normal_dnp_bam
            tumor_mapq0_bam=$tumor_mapq0_dnp_bam
            tumor_mapqnonzero_bam=$tumor_mapqnonzero_dnp_bam
        else
            normal_bam=$normal_snp_bam
            tumor_mapq0_bam=$tumor_mapq0_snp_bam
            tumor_mapqnonzero_bam=$tumor_mapqnonzero_snp_bam
        fi
        
        echo "Running mpileup on allele: $chr pos: $pos ..."
        samtools mpileup -r $region -f $fastafile -aa $normal_bam >> $normal_mpileupout
        samtools mpileup -r $region -f $fastafile -aa $tumor_mapq0_bam >> $tumor_mapq0_mpileupout
        samtools mpileup -r $region -f $fastafile -aa $tumor_mapqnonzero_bam >> $tumor_mapqnonzero_mpileupout
    fi

done < $somatic_mutations

awk '{gsub(/\$/, "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' $normal_mpileupout > ${normal_mpileupout/.txt/.tmp.txt}
mv ${normal_mpileupout/.txt/.tmp.txt} $normal_mpileupout
awk '{gsub(/\$/, "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' $tumor_mapq0_mpileupout > ${tumor_mapq0_mpileupout/.txt/.tmp.txt}
mv ${tumor_mapq0_mpileupout/.txt/.tmp.txt} $tumor_mapq0_mpileupout
awk '{gsub(/\$/, "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' $tumor_mapqnonzero_mpileupout > ${tumor_mapqnonzero_mpileupout/.txt/.tmp.txt}
mv ${tumor_mapqnonzero_mpileupout/.txt/.tmp.txt} $tumor_mapqnonzero_mpileupout


echo "Done!"


