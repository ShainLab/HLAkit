#!/bin/bash

set -euo pipefail

# Change MAPQ of multimapping reads to 0

usage() {
    cat << EOF
Usage: $0 --resultdir <dir> --alleles <path> --threads <int>

Required arguments:
    -o, --resultdir <dir>    Output directory for results
    -a, --alleles <file>     Updated alleles

Optional arguements:
    -t, --threads <int>      Number of threads
    -h, --help               Show this help message


Example:
    $0 --resultdir /path/to/resultdirectory --alleles allelelist_updated.txt
    $0 -o /path/to/resultdirectory -a allelelist_updated.txt
EOF
    exit 1
}

# Initialize variables
resultdir=""
alleles=""
threads=1

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--resultdir)
            resultdir="$2"
            shift 2
            ;;
        -a|--alleles)
            alleles="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
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
if [ -z "$resultdir" ] || [ -z "$alleles" ] ; then
    echo "Error: All arguments are required"
    usage
fi

# Validate resultdir and alleles file
if [ ! -d "$resultdir" ]; then
    echo "Error: Result directory does not exist: $resultdir"
    exit 1
fi
if [ ! -f "$alleles" ]; then
    echo "Error: Alleles file not found: $alleles"
    exit 1
elif [ ! -s "$alleles" ]; then
    echo "Error: Alleles file is empty: $alleles"
    exit 1
fi

if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: --threads must be a positive integer (got: '$threads')."
    exit 1
fi

RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m'
log_warn() {
    local msg="[$(date +'%Y-%m-%d %H:%M:%S')] [WARNING] $*"
    echo -e "${YELLOW}${msg}${NC}" >&2
    # echo "$msg" >> "$log_file"
}



echo "Change MAPQ for any read that maps to multiple alleles to 0 ..."
echo "Extracting readnames from:"
while read allele
do
    if [ -z "$allele" ]; then continue; fi
    for file in $resultdir/${allele}.*.[SD]NP.txt
    do
        if [ ! -f "$file" ]; then
            log_warn "Expected alignment file not found, skipping: $file"
            continue
        fi
        echo $file
        awk 'NR>4{print $1}' "$file" | sort | uniq  | egrep -v "^@" > "${file/.txt/.readnames.txt}"
    done
done < "$alleles"


for sampletype in tumor normal
do
    for muttype in SNP DNP
    do
        readname_files=("$resultdir"/hla_*.${sampletype}.${muttype}.readnames.txt)
        if [ "${#readname_files[@]}" -lt 2 ]; then
            log_warn "Warning: Only one readnames file found for $sampletype $muttype. Skipping readname comparison step. Creating empty bam files for MAPQ0 and copying ${readname_files/.readnames.txt/.txt} as the MAPQ corrected and non-zero bam file."
            cp ${readname_files/.readnames.txt/.bam} ${readname_files/.readnames.txt/.MAPQcorrected.bam}
            cp ${readname_files/.readnames.txt/.bam} ${readname_files/.readnames.txt/.MAPQnonzero.bam}
            samtools view -@ $threads -bS -H ${readname_files/.readnames.txt/.sam} > ${readname_files/.readnames.txt/.MAPQzero.bam}
            samtools index -@ $threads ${readname_files/.readnames.txt/.MAPQcorrected.bam}
            samtools index -@ $threads ${readname_files/.readnames.txt/.MAPQnonzero.bam}
            samtools index -@ $threads ${readname_files/.readnames.txt/.MAPQzero.bam}
        else
            while read allele; do
                if [ -z "$allele" ]; then continue; fi
                file="$resultdir"/${allele}.${sampletype}.${muttype}.txt

                if [ ! -f "$file" ]; then
                    log_warn "Alignment file not found for allele '$allele' ($sampletype $muttype): $file — skipping MAPQ correction for this allele."
                    continue
                fi

                echo "Currently working on: $file"

                echo "Extracting multimapping reads ... "
                tmpfile="$sampletype"."$muttype"."$allele".tmp.txt
                tmpmatchfile="$sampletype"."$muttype"."$allele".tmpmatch.txt
                cat $(ls "$resultdir"/hla_*.${sampletype}.${muttype}.readnames.txt | grep -v "^${file/.txt/.readnames.txt}$") > "$resultdir"/$tmpfile

                grep -Ff "$resultdir"/$tmpfile "$file" > "$resultdir"/$tmpmatchfile || true

                echo "Correcting MAPQ of multimapping reads ..."
                samtools view -@ $threads -H "${file/.txt/.sam}" > "${file/.txt/.MAPQcorrected.sam}"

                if [[ ! -s "$resultdir"/$tmpmatchfile ]]; then
                    echo "No multimapping reads found."
                    cp ${file/.txt/.bam} ${file/.txt/.MAPQcorrected.bam}
                    cp ${file/.txt/.bam} ${file/.txt/.MAPQnonzero.bam}
                    samtools view -@ $threads -H ${file/.txt/.sam} > ${file/.txt/.MAPQzero.sam}
                    samtools view -@ $threads -bS -o ${file/.txt/.MAPQzero.bam} ${file/.txt/.MAPQzero.sam} -> Changed this because the previous one was -H when the other code mentions -bS
                    samtools index -@ $threads ${file/.txt/.MAPQcorrected.bam}
                    samtools index -@ $threads ${file/.txt/.MAPQnonzero.bam}
                    samtools index -@ $threads ${file/.txt/.MAPQzero.bam}
                else
                    awk 'BEGIN{ FS = "\t"; OFS="\t" } NR==FNR { matches[$1]; next }
                         $1 in matches {
                             $5 = 0
                         }
                         { print }' "$resultdir"/$tmpmatchfile "${file}" >> "${file/.txt/.MAPQcorrected.sam}"

                    samtools sort -@ $threads -o "${file/.txt/.MAPQcorrected.sorted.sam}" "${file/.txt/.MAPQcorrected.sam}"
                    mv "${file/.txt/.MAPQcorrected.sorted.sam}" "${file/.txt/.MAPQcorrected.sam}"
                    samtools view -@ $threads -bS "${file/.txt/.MAPQcorrected.sam}" > "${file/.txt/.MAPQcorrected.bam}"
                    samtools index -@ $threads "${file/.txt/.MAPQcorrected.bam}"

                    echo "Creating separate bam files for MAPQzero and MAPQnonzero reads ..."
                    samtools view -@ $threads -H ${file/.txt/.sam} > "${file/.txt/.MAPQzero.sam}"
                    grep -v ^@ "${file/.txt/.MAPQcorrected.sam}" | awk '$5==0{print $0}' >> "${file/.txt/.MAPQzero.sam}"
                    samtools view -@ $threads -bS -o "${file/.txt/.MAPQzero.bam}" "${file/.txt/.MAPQzero.sam}"
                    samtools index -@ $threads "${file/.txt/.MAPQzero.bam}"

                    samtools view -@ $threads -H ${file/.txt/.sam} > "${file/.txt/.MAPQnonzero.sam}"
                    grep -v ^@ "${file/.txt/.MAPQcorrected.sam}" | awk '$5>0{print $0}' >> "${file/.txt/.MAPQnonzero.sam}"
                    samtools view -@ $threads -bS -o "${file/.txt/.MAPQnonzero.bam}" "${file/.txt/.MAPQnonzero.sam}"
                    samtools index -@ $threads "${file/.txt/.MAPQnonzero.bam}"
                fi
                echo "$file Done!"
            done < "$alleles"
        fi
    done
done

