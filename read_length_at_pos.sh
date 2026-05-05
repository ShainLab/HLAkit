#!/usr/bin/env bash

mml=$1
resultdir=$2

:>"$read_lengths_count"

while IFS='\t' read LINE
do
    read -r allele pos ref mut variant_classification <<< "$LINE"
    mut=${mut//$'\r'/}
    mut=$(echo "$mut" | tr -dc 'ACGTacgt')
    mut="${mut:0:1}"

    bam_file=$resultdir/novoalign_tumor.coordsort.dedup.RG.typedsubset.clean${variant_classification}.${allele}.sorted.bam

    read_lengths=$resultdir/tumor.${variant_classification}.${allele}.${ref}.${mut}.readlengths.txt
    read_lengths_count=$resultdir/tumor.${variant_classification}.${allele}.${ref}.${mut}.readlength_counts.txt

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
            print length($10)
        } else if (mut_type == "INS" && ("+" mut) == "+" read_base) {
            print length($10)
        } else if (mut_type == "DEL" && mut == "-") {
            print length($10)
        }
    }' > "$read_lengths"

    # Count read lengths
    fl=$(( $(samtools view $bam_file | awk '{print length($10)}' | sort -nu | tail -1) - 5 ))
    count_fulllength=$(awk -v fl=$fl '$1 >= fl { count++ } END { print count+0 }' "$read_lengths")
    count_short=$(awk -v fl=$fl '$1 < fl { count++ } END { print count+0 }' "$read_lengths")

    # Write final output:  allele, pos, mut, 101_count, <101_count, total
    printf "${allele}\t${pos}\t${ref}\t${mut}\t${count_fulllength}\t${count_short}\n" >> "$read_lengths_count"


done < "$mml"

