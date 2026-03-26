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

# Validate input has at least 2 columns (allele and POS)
ncols=$(awk 'NR==2{print NF; exit}' "$input")
if [ -z "$ncols" ] || [ "$ncols" -lt 2 ]; then
    echo "Error: Input file appears to have fewer than 2 columns (expected at least allele and POS). Got: ${ncols:-0} column(s)."
    exit 1
fi

# Validate bedfile has at least 4 columns (chrom, start, end, feature)
bed_ncols=$(awk 'NR==1{print NF; exit}' "$bedfile")
if [ -z "$bed_ncols" ] || [ "$bed_ncols" -lt 4 ]; then
    echo "Error: BED file appears to have fewer than 4 columns (expected chrom, start, end, feature). Got: ${bed_ncols:-0} column(s)."
    exit 1
fi

outfile="${input/.txt/.features.txt}"
awk 'NR==1{print $0"\tFeature"}' "$input" > "$outfile"
echo "Adding exon/intron/utr features to $input ..."

total=0
matched=0
while IFS= read -r LINE
do
	allele=$(echo $LINE | awk '{print $1}')
	POS=$(echo $LINE | awk '{print $2}')
    if [ "$allele" != "CHROM" ]; then
    	result=$(awk -v pos="$POS" -v line="$LINE" -v allele="$allele" 'NR>1 && $1==allele && $2<=pos && $3>=pos {print line"\t"$4}' $bedfile)
    	if [ -n "$result" ]; then
    		echo "$result" >> "$outfile"
    		matched=$((matched + 1))
    	else
    		echo "Warning: No feature found in BED file for allele '$allele' at POS $POS — row will be excluded from output."
    	fi
    	total=$((total + 1))
    fi
done <"$input"

if [ "$matched" -eq 0 ] && [ "$total" -gt 0 ]; then
    echo "Error: No mutations matched any feature in the BED file. Output file will contain only the header. Check that allele names and coordinates are consistent between input and BED file."
    exit 1
fi

echo "Features added successfully! ($matched/$total mutations annotated)"
