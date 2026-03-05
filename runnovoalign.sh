#!/bin/bash

set -euo pipefail

usage() {
    cat << EOF
Usage: $0 --f1 <path> --f2 <path> --threads <int> --nixFile <nix> --output <file> --hlakit <path> --format <string>
Required arguments:
    --f1 <path>            Path to subsetted normal or tumor fastq R1 file
    --f2 <path>            Path to subsetted normal or tumor fastq R2 file
    --threads <int>        Number of threads for parallel processing (should be a multiple of 4)
    --nixFile <nix>        Novoalign indexed reference file
    --output <file>        Name of output SAM file
    --hlakit <path>        Path to hlakit software

Optional arguements:
    --softClip <int>       Allow soft-clipping of reads. 0- do not allow; 1- allow (default: 0)
    --format <string>      FASTQ format (ILM1.8/ILMFQ/STDFQ; auto-detected if not specified)
    --help                 Show this help message


Example:
    $0 --f1 normal_subsetted.1.fastq --f2 normal_subsetted.2.fastq --threads 8 --nixFile hg19.nix --output novoalign_normal_reads.sam --hlakit $HOME/hlakit
EOF
    exit 1
}



# Initialize variables
f1=""
f2=""
threads=""
nixFile=""
output=""
hlakit=""
format=""
softClip=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --f1)
            f1="$2"
            shift 2
            ;;
        --f2)
            f2="$2"
            shift 2
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        --nixFile)
            nixFile="$2"
            shift 2
            ;;
        --output)
            output="$2"
            shift 2
            ;;
        --hlakit)
            hlakit="$2"
            shift 2
            ;;
        --format)
            format="$2"
            shift 2
            ;;
        --softClip)
            softClip="$2"
            shift 2
            ;;
        --help)
            usage
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [ -z "$f1" ] || [ -z "$f2" ] || [ -z "$threads" ] ||[ -z "$nixFile" ] || [ -z "$output" ] || [ -z "$hlakit" ]; then
    echo "Error: All arguments are required"
    usage
fi

# Validate files exist
for file in "$f1" "$f2" "$nixFile"; do
    if [ ! -f "$file" ]; then
        echo "Error: File file not found: $file"
        exit 1
    elif [ ! -s "$file" ]; then
        echo "Error: File is empty: $file"
        exit 1
    fi
done

echo "Running: runnovoalign.sh --f1 $f1 --f2 $f2 --threads $threads --nixFile $nixFile --output $output --hlakit $hlakit --format $format --softClip $softClip"


digits=2


novoalignPath=$hlakit/binaries/novoalign

if [ -z $format ]; then
    set +e
    echo "Fastq format not provided in input. Evaluating format of fastq files ..."
    qual_scores=$(awk 'NR % 4 == 0' "$f1" | head -n 10000)
    
    # Get ASCII values of quality characters
    min_qual=$(echo "$qual_scores" | \
        od -An -tuC | \
        awk '{for(i=1;i<=NF;i++) if($i>32) print $i}' | \
        sort -n | head -1)

    max_qual=$(echo "$qual_scores" | \
        od -An -tuC | \
        awk '{for(i=1;i<=NF;i++) if($i>32) print $i}' | \
        sort -n | tail -1)
    # Illumina 1.8+ (Phred+33): ASCII 33-73 (quality 0-40)
    # Sanger/Standard (Phred+33): ASCII 33-126 (quality 0-93)
    
    if [ "$min_qual" -ge 33 ] && [ "$max_qual" -le 74 ]; then
        format=ILM1.8
    elif [ "$min_qual" -ge 33 ] && [ "$max_qual" -le 126 ]; then
        format=STDFQ
    elif [ "$min_qual" -ge 64 ]; then
        format=ILMFQ
    else
        format=STDFQ
    fi

    echo "Fastq format determined as: $format."
fi

set -e

if [ $threads -gt 1 ]
    then echo "Running novoalign with $threads parallel jobs..."
else
    echo "Running novoalign using single thread ..."
fi

# Get line count
n1=$(wc -l "$f1" | awk '{print $1}')
n2=$(wc -l "$f2" | awk '{print $1}')

if [ "$n1" -ne "$n2" ]; then
    echo "Number of reads are not equal in R1 and R2 fastq files." >&2
    exit 1
fi

echo "Calculating number of reads for each split fastq file ..."
nPerFile=$(awk -v n1="$n1" -v threads="$threads" 'BEGIN {
    nPerFile = n1/threads;
    ceil_val = int(nPerFile);
    if (ceil_val != nPerFile) {
        nPerFile = int(nPerFile) + 1;
    }
    while (nPerFile % 4 != 0) {
        nPerFile = nPerFile + 1;
    }
    print nPerFile;
}')

ceil=$(awk -v n1="$n1" -v threads="$threads" 'BEGIN {
    nPerFile = n1/threads;
    print int(nPerFile) + (nPerFile > int(nPerFile) ? 1 : 0);
}')

echo "Total reads: F1=$n1 F2=$n2"
echo "Reads per split fastq file= $nPerFile"

split -l "$nPerFile" -d -a 2 "$f1" "$f1"
split -l "$nPerFile" -d -a 2 "$f2" "$f2"

echo "Aligning the split fastq files ..."
for ((i = 0; i < threads; i++)); do
    (
        # Format suffix with leading zeros
        suffix=$(printf "%0${digits}d" "$i")
        
        
        if [ "$softClip" -eq 0 ]; then
            "$novoalignPath" -d "$nixFile" -f "${f1}${suffix}" "${f2}${suffix}" \
                -F "$format" -R 0 -r all -o SAM -o FullNW | grep -P '\thla' > "${output}${suffix}"
        else
            "$novoalignPath" -d "$nixFile" -f "${f1}${suffix}" "${f2}${suffix}" \
                -F "$format" -R 0 -r all -o SAM -g 20 -x 3 | grep -P '\thla' > "${output}${suffix}"
        fi
    ) &
done

wait

echo "Concatenating SAM files into one ..."
for ((i = 0; i < threads; i++)); do
    suffix=$(printf "%0${digits}d" "$i")
    cat "${output}${suffix}" >> "$output"
    rm -rf "${output}${suffix}" "${f1}${suffix}" "${f2}${suffix}"
done

echo "Novoalign run finished. Output written to: $output."

