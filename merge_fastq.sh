#!/bin/bash

set -euo pipefail

usage() {
    cat << EOF
Usage: $0 --kmer <file> --chr6region <file> --output <file>

Required arguments:
    -k, --kmer <file>         kmer subsetted fastq file
    -c, --chr6region <file>   chr6 HLA region subsetted fastq file
    -o, --output <file>       output file

Optional arguments:
    -h, --help                Show this help message


Example:
    $0 --kmer tumor_kmer --chr6region chr6region.tumor --output tumor_subsetted
    $0 -k tumor_kmer -c chr6region.tumor -o tumor_subsetted
EOF
    exit 1
}

# Initialize variables
kmer=""
chr6region=""
output=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -k|--kmer)
            kmer="$2"
            shift 2
            ;;
        -c|--chr6region)
            chr6region="$2"
            shift 2
            ;;
        -o|--output)
            output="$2"
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
if [ -z "$kmer" ] || [ -z "$chr6region" ] || [ -z "$output" ]; then
    echo "Error: All arguments are required"
    usage
fi

# Validate files exist
for file in "$kmer".1.fastq "$kmer".2.fastq "$chr6region".1.fastq "$kmer".2.fastq ; do
    if [ ! -f "$file" ]; then
        echo "Error: File not found: $file"
        exit 1
    elif [ ! -s "$file" ]; then
        echo "WARNING: File is empty: $file"
    fi
done


echo "Merging R1 reads ..."
cat "${kmer}.1.fastq" "${chr6region}.1.fastq" | \
awk 'NR%4==1 {h=$0; getline s; getline p; getline q
     idx=index(h,"/"); name=substr(h,1,idx-1)
     if(!seen[name]++) print h"\n"s"\n"p"\n"q
}' > "${output}.1.fastq"
echo "Merging R1 reads done."

echo "Merging R2 reads ..."
cat "${kmer}.2.fastq" "${chr6region}.2.fastq" | \
awk 'NR%4==1 {h=$0; getline s; getline p; getline q
     idx=index(h,"/"); name=substr(h,1,idx-1)
     if(!seen[name]++) print h"\n"s"\n"p"\n"q
}' > "${output}.2.fastq"

echo "Merging R2 reads done."

echo "Done! Merged $(wc -l < ${output}.1.fastq | awk '{print $1/4}') unique reads"


