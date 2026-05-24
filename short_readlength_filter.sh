#!/usr/bin/env bash

usage() {
    cat << EOF
Usage: $0 --resultdir <dir>

Required arguments:
    -o, --resultdir <dir>  Output directory for results

Example:
    $0 --resultdir /path/to/resultdirectory
EOF
    exit 1
}

# Initialize variables
resultdir=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--resultdir)
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

if [ -z "$resultdir" ]; then
    echo "Error: All arguments are required"
    usage
fi


mkdir -p $resultdir/shortread_filter

shortreadlength_filter=$resultdir/shortread_filter/shortreadlength_filter.txt
:> "$shortreadlength_filter"
 
echo "Running short read length filter ..."
 

while IFS=$'\t' read LINE
do
    # read -r CHROM POS REF ALT FILTER Normal_Ref Normal_Mut Tumor_Ref Tumor_Mut Tumor_MAF UV Feature amino_acid_change variant_classification CDS_pos aa_pos variant_type Artifacts Zygosity tumor_purity zcMAF  <<< "$LINE"
    CHROM=$(echo "$LINE" | cut -f1)
    POS=$(echo "$LINE" | cut -f2)
    REF=$(echo "$LINE" | cut -f3)
    ALT=$(echo "$LINE" | cut -f4)
    variant_type=$(echo "$LINE" | cut -f17)

    ALT=${ALT//$'\r'/}
    ALT=$(echo "$ALT" | tr -dc 'ACGTacgt')
    variant_type=$(echo "$variant_type" | tr -dc 'SNPDNPINSDEL')
    
    if [[ "$variant_type" == "DNP" ]]
    then
        ALT="${ALT:0:1}"
        bam_file=$resultdir/${CHROM}.tumor.DNP.bam
    else
        bam_file=$resultdir/${CHROM}.tumor.SNP.bam
    fi

    if [ -z "$ALT" ]; then
        echo "Warning: Could not parse mutation for ${CHROM}:${POS} ${REF}>${ALT} — skipping."
        printf $CHROM"\t"$POS"\t"$REF"\t"$ALT"\t\n"  >> $shortreadlength_filter
        continue
    fi

    if [ ! -f "$bam_file" ]; then
        echo "Warning: BAM file not found: $bam_file — skipping ${CHROM}:${POS} ${REF}>${ALT}"
        printf $CHROM"\t"$POS"\t"$REF"\t"$ALT"\t\n"  >> $shortreadlength_filter
        continue
    fi

    if [ ! -f "${bam_file}.bai" ] && [ ! -f "${bam_file%.bam}.bai" ]; then
        echo "Warning: No BAM index found for $bam_file — samtools view may fail."
    fi

    echo "Current mutation: "
    printf $CHROM"\t"$POS"\t"$REF"\t"$ALT"\t\n" 

    # Extract reads overlapping the position using samtools
    output_file=$resultdir/shortread_filter/${CHROM}.${POS}.${ALT}.${variant_type}.readlengths.txt
    if [[ "$variant_type" == "DEL" ]] || [[ "$variant_type" == "INS" ]]; then
        ALT_TAG=${ALT:0:10}
        output_file=$resultdir/shortread_filter/${CHROM}.${POS}.${ALT_TAG}.${variant_type}.readlengths.txt
    fi
    
    if [[ "$variant_type" == "DNP" ]] || [[ "$variant_type" == "SNP" ]]; then
        :> "$output_file"
        echo "Checking read lengths in mutation: ${CHROM}:${POS} ${REF}>${ALT}"
        samtools view "$bam_file" "$CHROM:$POS-$POS" | \
        awk -v pos="$POS" -v mut="$ALT" '
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
        }' > "$output_file"
    fi

    # Count read lengths
    full_length=$(samtools view $bam_file | head -500 | awk '{print length($10)}' | sort -n | tail -1)
    
    if [ -z "$full_length" ] || [ "$full_length" -eq 0 ]; then
        echo "Warning: Could not determine full read length from $bam_file — skipping ${CHROM}:${POS} ${REF}>${ALT}"
        printf $CHROM"\t"$POS"\t"$REF"\t"$ALT"\t\n"  >> $shortreadlength_filter
        continue
    fi

    count_full=$(awk -v full_length=$full_length '$1 >= full_length { count++ } END { print count+0 }' "$output_file")
    count_short=$(awk -v full_length=$full_length '$1 < full_length { count++ } END { print count+0 }' "$output_file")

    echo "Full length reads in ${CHROM}:${POS} ${REF}>${ALT} =" $count_full
    echo "Short length reads in ${CHROM}:${POS} ${REF}>${ALT} =" $count_short

    if [[ "$count_full" -eq 0 ]]; then
        printf $CHROM"\t"$POS"\t"$REF"\t"$ALT"\tartifact_shortreadlength\n" >> $shortreadlength_filter
    else
        printf $CHROM"\t"$POS"\t"$REF"\t"$ALT"\t\n"  >> $shortreadlength_filter
    fi
    
done < <(tail -n +2 $resultdir/hla_somatic_mutations.txt)

if [ ! -s "$shortreadlength_filter" ]; then
    echo "Warning: shortreadlength_filter is empty — no mutations were processed."
fi

awk -F"\t" 'BEGIN { OFS="\t" }
FNR==NR { key=$1 FS $2 FS $3; f2_col5[key]=$5; next }
FNR==1 { print; next }
{
    key=$1 FS $2 FS $3
    if (key in f2_col5 && f2_col5[key]=="artifact_shortreadlength") {
        if ($18=="") $18="artifact_shortreadlength"
        else $18=$18",artifact_shortreadlength"
    }
    print
}' "$shortreadlength_filter" "$resultdir/hla_somatic_mutations.txt" > "$resultdir/hla_somatic_mutations.temp.txt"

# rm -r $resultdir/shortread_filter

if [ ! -s "$resultdir/hla_somatic_mutations.temp.txt" ]; then
    echo "Error: Output mutation file is empty after short read filter merge — not overwriting original."
else
    mv "$resultdir/hla_somatic_mutations.temp.txt" "$resultdir/hla_somatic_mutations.txt"
fi

echo "Short read length filter done."


