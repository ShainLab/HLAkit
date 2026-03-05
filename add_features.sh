#!/bin/bash

set -euo pipefail

# add feature information (exon, intron, utr) to somatic mutation files

usage() {
    cat << EOF
Usage: $0 --input <file> --bedfile <file>

Required arguments:
    -i, --input <file>      Input somatic mutation file
    -b, --bedfile <file>     bedfile with feature information

Optional arguements:
    -h, --help                Show this help message


Example:
    $0 --input somatic_mutations.UV.txt --bedfile hla_bedfile.bed
    $0 -i somatic_mutations.UV.txt -b hla_bedfile.bed
EOF
    exit 1
}


# Initialize variables
input=""
bedfile=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            input="$2"
            shift 2
            ;;
        -b|--bedfile)
            bedfile="$2"
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
if [ -z "$input" ] || [ -z "$bedfile" ]; then
    echo "Error: All arguments are required"
    usage
fi

# Validate files exist
for file in $input $bedfile
do
	if [ ! -f "$file" ]; then
	    echo "Error: File not found: $file"
	    exit 1
	elif [ ! -s "$file" ]; then
	    echo "Error: File is empty: $file"
	    exit 1
	fi
done


awk 'NR==1{print $0"\tFeature"}' "$input" > "${input/.txt/.features.txt}"

echo "Adding exon/intron/utr features to $input ..."
while IFS= read -r LINE
do
	allele=$(echo $LINE | awk '{print $1}')
	POS=$(echo $LINE | awk '{print $2}')

	awk -v pos="$POS" -v line="$LINE" -v allele="$allele" 'NR>1 && $1==allele && $2<=pos && $3>=pos {print line"\t"$4}' $bedfile >> ${input/.txt/.features.txt}

done <"$input"

echo "Features added successfully!"

