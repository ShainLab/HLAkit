#!/bin/bash

set -euo pipefail


usage() {
    cat << EOF
Usage: $0 --resultdir <path> --allelelist <file> --hlakit <path>

Required arguments:
    -o, --resultdir <path>    Path to fastq file
    -a, --allelelist <file>   Updated allele list file
    -p, --hlakit <path>  Path to hlakit repo

Optional arguements:
    -h, --help                Show this help message


Example:
    $0 --resultdir resultdir --allelelist allelelist.updated.txt --hlakit /path/to/hlakit
EOF
    exit 1
}

# Initialize variables
resultdir=""
allelelist=""
hlakit=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--resultdir)
            resultdir="$2"
            shift 2
            ;;
        -a|--allelelist)
            allelelist="$2"
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
if [ -z "$resultdir" ] || [ -z "$allelelist" ] || [ -z "$hlakit" ]; then
    echo "Error: All arguments are required"
    usage
fi

# Validate files exist
if [ ! -f "$allelelist" ]; then
    echo "Error: File not found: $allelelist"
    exit 1
elif [ ! -s "$allelelist" ]; then
    echo "Error: File is empty: $allelelist"
    exit 1
fi


# "Run hs metrics for allele specific HLA bam file"
interval_list=$hlakit/resources/hs_metrics_interval_list.txt
if [ ! -f "$interval_list" ]; then
    echo "Error: File not found: $interval_list"
    exit 1
fi

while read -r allele
do
    allele=$(echo "$allele" | tr -d '\r' | xargs)
    echo "Creating interval file for allele: $allele ..."
    baits=$resultdir/"$allele".interval_list.txt
    egrep "\b${allele}\b" "$interval_list" > "$baits"


    bamfile_list=$(ls $resultdir/${allele}.*.SNP.MAPQ*zero.bam)

    if [ -n "$bamfile_list" ]
    then
        for bamfile in $resultdir/${allele}.*.SNP.MAPQ*zero.bam
        do
            echo "Running HS metrics on "$bamfile" ..."
            gatk CollectHsMetrics \
                -BI "$baits" \
                -TI "$baits" \
                -I "$bamfile" \
                -O ${bamfile/.bam/.hsmetrics.txt}
        done

        echo "Calculating coverage in MAPQ zero and nonzero bam files ..."
        mapq0_normal=$(awk 'NR==8{print $10}' $resultdir/${allele}.normal.SNP.MAPQzero.hsmetrics.txt)
        mapqnonzero_normal=$(awk 'NR==8{print $10}' $resultdir/${allele}.normal.SNP.MAPQnonzero.hsmetrics.txt)
        mapq0_tumor=$(awk 'NR==8{print $10}' $resultdir/${allele}.tumor.SNP.MAPQzero.hsmetrics.txt)
        mapqnonzero_tumor=$(awk 'NR==8{print $10}' $resultdir/${allele}.tumor.SNP.MAPQnonzero.hsmetrics.txt)

        :> $resultdir/"$allele".coverage.txt
        normal_cov=$(awk -v n0="${mapq0_normal:-0}" -v n1="${mapqnonzero_normal:-0}" 'BEGIN {print (n0 / 2) + n1}')
        tumor_cov=$(awk -v t0="${mapq0_tumor:-0}" -v t1="${mapqnonzero_tumor:-0}" 'BEGIN {print (t0 / 2) + t1}')
        printf "$allele.normal\t""$normal_cov\n" >> $resultdir/"$allele".coverage.txt
        printf "$allele.tumor\t""$tumor_cov\n" >> $resultdir/"$allele".coverage.txt
    fi
done < "$allelelist"

echo "Done!"

