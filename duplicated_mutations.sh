#!/bin/bash

# Find duplicate mutations in 2 alleles


usage() {
    cat << EOF
Usage: $0 --mml <file>

Required arguments:
    -m, --mml <file>          Path to mml file
    -r, --resultdir <path>    Result directory

Optional arguements:
    -h, --help                Show this help message


Example:
    $0 --mml mml
EOF
    exit 1
}

# Initialize variables
mml=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--mml)
            mml="$2"
            shift 2
            ;;
        -r|--resultdir)
            resultdir="$2"
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
if [ -z "$mml" ] || [ -z "$resultdir" ] ; then
    echo "Error: All arguments are required"
    usage
fi

# Validate files exist
if [ ! -f "$mml" ]; then
    echo "Error: File not found: $mml"
    exit 1
elif [ ! -s "$mml" ]; then
    echo "Error: File is empty: $mml"
    exit 1
fi

# Validate resultdir exists
if [ ! -d "$resultdir" ]; then
    echo "Error: Result directory does not exist: $resultdir"
    exit 1
fi

# Validate samtools is available
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found in PATH."
    exit 1
fi

threshold=4

echo "Looking for duplicate mutations ..."

awk -F"\t" '
{
    rows[NR] = $0
    f1[NR] = $1
    f1p[NR] = substr($1, 1, 5)
    f2[NR] = $2
    f4[NR] = $4
    f5[NR] = $5
    f19[NR] = $19
    total = NR
}
END {
    for (i = 1; i <= total; i++) {
        if (used[i]) continue
        matched = 0

        # first pass: exact col2 match
        for (j = i+1; j <= total; j++) {
            if (used[j]) continue
            if (f1[i] == f1[j]) continue
            if (f1p[i] != f1p[j]) continue
            if (f5[i] != f5[j]) continue
            if (f19[i] != f19[j]) continue
            if (length(f4[i]) != length(f4[j])) continue
            if (f2[i] != f2[j]) continue
            used[i] = 1
            used[j] = 1
            split(rows[i], r1)
            split(rows[j], r2)
            print r1[1] "\t" r2[1] "\t" r1[2] "\t" r2[2] "\t" r1[4] "\t" r2[4] "\t" r1[5] "\t" r1[24]
            matched = 1
            break
        }

        if (matched) continue

        # second pass: within 200 bases
        for (j = i+1; j <= total; j++) {
            if (used[j]) continue
            if (f1[i] == f1[j]) continue
            if (f1p[i] != f1p[j]) continue
            if (f5[i] != f5[j]) continue
            if (f19[i] != f19[j]) continue
            if (length(f4[i]) != length(f4[j])) continue
            pos_diff = f2[i] - f2[j]
            if (pos_diff < 0) pos_diff = -pos_diff
            if (pos_diff > 200) continue
            used[i] = 1
            used[j] = 1
            split(rows[i], r1)
            split(rows[j], r2)
            print r1[1] "\t" r2[1] "\t" r1[2] "\t" r2[2] "\t" r1[4] "\t" r2[4] "\t" r1[5] "\t" r1[24]
            break
        }
    }
}
' FS="\t" $mml > $resultdir/mutations_to_dedup.txt

if [ ! -s "$resultdir/mutations_to_dedup.txt" ]; then
    echo "No duplicate mutations found. Nothing to deduplicate. Exiting."
    cp $mml ${mml/.txt/.dedup.txt}
    exit 0
fi

echo "Found $(wc -l < $resultdir/mutations_to_dedup.txt) mutation pair(s) to deduplicate."

# :>$resultdir/duplicated_mutations.txt

cp $mml ${mml/.txt/.dedup.txt}
mml=${mml/.txt/.dedup.txt}

while read LINE
do
    read allele1 allele2 pos1 pos2 ref1 ref2 mut variant_type <<< "$LINE"
    mut=${mut//$'\r'/}
    mut=$(echo "$mut" | tr -dc 'ACGTacgt')
    mut="${mut:0:1}"

    if [ -z "$mut" ]; then
        echo "Warning: Could not parse a valid nucleotide for mutation in line: $LINE — skipping."
        continue
    fi

    if [ "$variant_type" != "DNP" ]; then variant_type=SNP; fi

    for num in 1 2
    do
        allele="allele${num}" ; allele="${!allele}"
        pos="pos${num}" ; pos="${!pos}"
        ref="ref${num}" ; ref="${!ref}"

        bam_file=$resultdir/${allele}.tumor.${variant_type}.bam

        if [ ! -f "$bam_file" ]; then
            echo "Error: BAM file not found for allele ${allele} (${variant_type}): $bam_file"
            exit 1
        fi
        if [ ! -f "${bam_file}.bai" ] && [ ! -f "${bam_file%.bam}.bai" ]; then
            echo "Warning: No BAM index (.bai) found for $bam_file — samtools view may fail."
        fi

        readnames=$resultdir/${allele}.${pos}.${ref}.${mut}.${variant_type}.readnames_dup_mut_check.txt

        # Extract reads overlapping the position using samtools
        samtools view "$bam_file" "$allele:$pos-$pos" | \
        awk -v pos="$pos" -v mut="$mut" '
        function get_base_and_type(seq, cig, start, pos,  ref_pos, read_pos, len, op) {
            ref_pos = start
            read_pos = 1
            while (match(cig, /[0-9]+[MIDNSHP=X]/)) {
                len = substr(cig, RSTART, RLENGTH - 1)
                op = substr(cig, RSTART + length(len), 1)
                cig = substr(cig, RSTART + RLENGTH)
                if (op == "M" || op == "=" || op == "X") {
                    if (pos >= ref_pos && pos < ref_pos + len) {
                        base = substr(seq, read_pos + (pos - ref_pos), 1)
                        return base "|" "SNP"
                    }
                    ref_pos += len
                    read_pos += len
                } else if (op == "I") {
                    if (ref_pos == pos) {
                        ins = substr(seq, read_pos, len)
                        return ins "|" "INS"
                    }
                    read_pos += len
                } else if (op == "D") {
                    if (ref_pos == pos) {
                        return "-" "|" "DEL"
                    }
                    ref_pos += len
                } else if (op == "N") {
                    ref_pos += len
                } else if (op == "S") {
                    read_pos += len
                } else if (op == "H") {
                    continue
                }
            }
            return "" "|" "NONE"
        }
        {
            split(get_base_and_type($10, $6, $4, pos), result, "|")
            read_base = result[1]
            mut_type = result[2]

            if (mut_type == "SNP" && read_base == mut) {
                print $1
            } else if (mut_type == "INS" && ("+" mut) == "+" read_base) {
                print $1
            } else if (mut_type == "DEL" && mut == "-") {
                print $1
            }
        }' > "$readnames"
    done

    file1=$resultdir/${allele1}.${pos1}.${ref1}.${mut}.${variant_type}.readnames_dup_mut_check.txt
    file2=$resultdir/${allele2}.${pos2}.${ref2}.${mut}.${variant_type}.readnames_dup_mut_check.txt

    if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
        echo "Error: Read name file(s) missing after samtools step — file1: $file1, file2: $file2. Exiting."
        exit 1
    fi

    total1=$(wc -l < $file1)
    total2=$(wc -l < $file2)
    common_reads=$(comm -12 <(sort $file1) <(sort $file2) | wc -l)
    echo "Total mut reads in ${allele1}:${pos1}:${ref1}>${mut} = $total1"
    echo "Total mut reads in ${allele2}:${pos2}:${ref2}>${mut} = $total2"
    echo "Common reads = $common_reads"
    is_duplicate=0

    thresh=0.10
    if [ $total1 -gt 0 ] && [ $total2 -gt 0 ]; then
        pct1=$(echo "scale=4; $common_reads / $total1" | bc -l)
        pct2=$(echo "scale=4; $common_reads / $total2" | bc -l)
        if (( $(echo "$pct1 >= $thresh" | bc -l) )) || (( $(echo "$pct2 >= $thresh" | bc -l) )); then
            is_duplicate=1
        fi
    fi

    if [ $is_duplicate -eq 0 ]; then
        echo "${allele1}:${pos1}:${ref1}>${mut} and ${allele2}:${pos2}:${ref2}>${mut} are unique mutations."
    fi


    if [ $is_duplicate -eq 1 ]; then
        echo "Mutation with higher MAF will be reported. If both mutations have same MAF, HLAkit will return the first mutation."

        allele1MAF=$(awk -v a=$allele1 -v p=$pos1 -v r=$ref1 -v m=$mut '$1==a && $2==p && $4==r && $5==m{print $16}' $mml)
        allele2MAF=$(awk -v a=$allele2 -v p=$pos2 -v r=$ref2 -v m=$mut '$1==a && $2==p && $4==r && $5==m{print $16}' $mml)

        if [ -z "$allele1MAF" ] || [ -z "$allele2MAF" ]; then
            echo "Warning: Could not retrieve MAF for one or both alleles (allele1MAF='$allele1MAF', allele2MAF='$allele2MAF') — check that columns 1,2,4,5,16 in the mml file match expected format. Skipping deduplication for this pair."
            continue
        fi

        if [ "$allele1MAF" == "$allele2MAF" ]; then
            awk -v a=$allele2 -v p=$pos2 -v r=$ref2 -v m=$mut '!($1 == a && $2 == p && $4 == r && $5 == m)' $mml > ${mml/.txt/.tmp.txt}
        elif (( $(echo "$allele1MAF > $allele2MAF" | bc -l) )); then
            awk -v a=$allele2 -v p=$pos2 -v r=$ref2 -v m=$mut '!($1 == a && $2 == p && $4 == r && $5 == m)' $mml > ${mml/.txt/.tmp.txt}
        elif (( $(echo "$allele1MAF < $allele2MAF" | bc -l) )); then
            awk -v a=$allele1 -v p=$pos1 -v r=$ref1 -v m=$mut '!($1 == a && $2 == p && $4 == r && $5 == m)' $mml > ${mml/.txt/.tmp.txt}
        fi

        if [ ! -s "${mml/.txt/.tmp.txt}" ]; then
            echo "Error: Temp mml file is empty after deduplication awk step for ${allele1}:${pos1} / ${allele2}:${pos2} — awk filter may have removed all rows. Not overwriting mml."
            exit 1
        fi

        mv ${mml/.txt/.tmp.txt} $mml
    fi

done < $resultdir/mutations_to_dedup.txt

echo "Mutation deduplication done!"
