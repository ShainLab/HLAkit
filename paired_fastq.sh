#!/bin/bash

set -euo pipefail


usage() {
    cat << EOF
Usage: $0 --fastqfile <file>

Required arguments:
    -f, --fastqfile <file>    Path to fastq file

Optional arguements:
    -h, --help                Show this help message


Example:
    $0 --fastqfile fastqfile
EOF
    exit 1
}

# Initialize variables
fastqfile=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--fastqfile)
            fastqfile="$2"
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
if [ -z "$fastqfile" ]; then
    echo "Error: All arguments are required"
    usage
fi

# Validate files exist
if [ ! -f "$fastqfile" ]; then
    echo "Error: File not found: $fastqfile"
    exit 1
elif [ ! -s "$fastqfile" ]; then
    echo "Error: File is empty: $fastqfile"
    exit 1
fi


echo "Removing unpaired fastq reads from $fastqfile file ..."

outfile="${fastqfile}.paired"

# read 4 lines at a time
awk 'NR%4==1 && /\// {print; for(i=1;i<=3;i++) {getline; print}}' "$fastqfile" > "$outfile"


mv "$outfile" "$fastqfile"
echo "Done!"
