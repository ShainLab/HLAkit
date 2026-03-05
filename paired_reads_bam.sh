#!/bin/bash

set -euo pipefail


usage() {
    cat << EOF
Usage: $0 --input <file>

Required arguments:
    -i, --input <file>    Path to bam file

Optional arguements:
    -h, --help            Show this help message


Example:
    $0 --input bamfile.bam
EOF
    exit 1
}

# Initialize variables
input=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            input="$2"
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
if [ -z "$input" ]; then
    echo "Error: All arguments are required"
    usage
fi

# Validate files exist
if [ ! -f "$input" ]; then
    echo "Error: File not found: $input"
    exit 1
elif [ ! -s "$input" ]; then
    echo "Error: File is empty: $input"
    exit 1
fi



awk '
{
    if (buffered) {
        line1 = buffer
        split(line1, f1, "\t")
        line2 = $0
        split($0, f2, "\t")
        buffered = 0
    } else {
        line1 = $0
        split($0, f1, "\t")
        if (getline line2 <= 0) exit
        split(line2, f2, "\t")
    }
    
    if (f1[1] "" == f2[1] "") {
        print line1
        print line2
    } else {
        buffer = line2
        buffered = 1
    }
}
' "$input" > "${input/.txt/.temp.txt}"
mv "${input/.txt/.temp.txt}" "$input"

