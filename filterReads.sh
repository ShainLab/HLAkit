#!/bin/bash
set -euo pipefail
# filter non primary alignment, set mapping quality to 70, and set mismatch threshold to mm
usage() {
    cat << EOF
Usage: $0 --input <file> --output <file> --mismatch <int>
Required arguments:
    -i, --input <file>      Input BAM file for filtering  
    -o, --output <file>     Filename for filtered output BAM
    -m, --mismatch <int>    Maximum number of mismatches allowed per read
Optional arguements:
    -h, --help                Show this help message
Example:
    $0 --input novoalign_normal.coordsort.dedup.RG.bam --output novoalign_normal.coordsort.dedup.RG.cleanSNP.bam --mismatch 1
    $0 -i novoalign_normal.coordsort.dedup.RG.bam -o novoalign_normal.coordsort.dedup.RG.cleanSNP.bam -m 1
EOF
    exit 1
}
# Initialize variables
input=""
output=""
mismatch=""
# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            input="$2"
            shift 2
            ;;
        -o|--output)
            output="$2"
            shift 2
            ;;
        -m|--mismatch)
            mismatch="$2"
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
if [ -z "$input" ] || [ -z "$output" ] || [ -z "$mismatch" ]; then
    echo "Error: All arguments are required"
    usage
fi

# Validate mismatch is a non-negative integer
if ! [[ "$mismatch" =~ ^[0-9]+$ ]]; then
    echo "Error: --mismatch must be a non-negative integer (got: '$mismatch')."
    exit 1
fi

# Validate files exist
if [ ! -f "$input" ]; then
    echo "Error: File not found: $input"
    exit 1
elif [ ! -s "$input" ]; then
    echo "Error: File is empty: $input"
    exit 1
fi

echo "Processing the BAM file to filter out reads with Not Primary Alignment and more than threshold ($mismatch) mismatches."
echo "Changing the MAPQ to 70 to run Mutect2 later."
cat "$input" | awk -v mismatch="$mismatch" '
BEGIN {
    FS=OFS="\t"
}
# Print headers as-is
/^@/ {
    print
    next
}
# Process alignment lines
{
    # Skip if mates dont map to same reference
    if ($7 != "=") next
    
    # Skip unmapped reads (FLAG bit 0x4) as they error gatk CollectHsMetrics
    if (and($2, 4)) next
    
    cigar = $6
    
    # Parse CIGAR to count insertions and deletions (bases, not events)
    insertionCount = 0
    deletionCount = 0
    n = split(cigar, cigArr, /[MIDNSHPX=]/)
    for (i = 1; i < n; i++) {
        op = substr(cigar, length(cigArr[i]) + sum_prev + 1, 1)
        if (op == "I") insertionCount += cigArr[i]
        if (op == "D") deletionCount += cigArr[i]
        sum_prev += length(cigArr[i]) + 1
    }
    sum_prev = 0
    
    # Extract NM tag (edit distance)
    editCount = 0
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^NM:i:/) {
            split($i, nm, ":")
            editCount = nm[3]
            break
        }
    }
    
    # Count events (occurrences, not bases)
    insertionEventCount = gsub(/I/, "I", cigar)
    deletionEventCount = gsub(/D/, "D", cigar)
    
    # Calculate mismatch events
    mismatchEventCount = editCount - insertionCount - deletionCount
    
    # Total events
    eventCount = insertionEventCount + deletionEventCount + mismatchEventCount
    
    # Filter by event mismatch
    if (eventCount > mismatch) next
    # if (eventCount != 0 && eventCount != mismatch) next
    
    # Clear bit 8 (256 = not primary alignment)
    $2 = and($2, compl(256))
    
    # Set mapping quality to 70
    $5 = 70
    
    print
}
' >  "$output"

if [ ! -s "$output" ]; then
    echo "Warning: Output file is empty after filtering: $output — all reads were filtered out. Check that --mismatch $mismatch is not too strict for this input, and that reads have NM tags and mate reference '=' in the input file."
fi

echo "Finished filtering the BAM file."
echo "Output: $output"

