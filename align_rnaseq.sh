#!/bin/bash
set -uo pipefail
# generate RNA-Seq bam files
usage() {
    cat << EOF
Usage: $0 --inputbam <file> --resultdir <string> --ref <file> --gtf <file> --kmer <file> --allelefile <file> --reference <hg19|hg38>
Required arguments:
    -i, --inputbam <file>        Input BAM file for filtering
    -o, --resultdir <string>     Output directory for results
    -r, --ref <file>             Reference fasta file for HLA alleles
    -g, --gtf <file>             Annotation file in GTF format
    -k, --kmer <file>            kmer file for HLA alleles
    -a, --allelefile <file>      Path to allele typing result file
    --reference <string>         GRCh reference build used to align the bam files. Accepted values are hg19 and hg38
Optional arguments:
    -t, --threads <int>          Number of threads to use (default: 1)
    -h, --help                   Show this help message
Example:
    $0 --inputbam rnaseq.bam --resultdir /path/to/resultdirectory --ref /path/to/hla.fasta --gtf /path/to/hla.gtf --kmer /path/to/kmer.txt --allelefile /path/to/allelefile --reference hg38
EOF
    exit 1
}

# Initialize variables
inputbam=""
resultdir=""
ref=""
gtf=""
kmer=""
allelefile=""
reference=""
threads=1

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--inputbam) inputbam="$2"; shift 2 ;;
        -o|--resultdir) resultdir="$2"; shift 2 ;;
        -r|--ref) ref="$2"; shift 2 ;;
        -g|--gtf) gtf="$2"; shift 2 ;;
        -k|--kmer) kmer="$2"; shift 2 ;;
        -a|--allelefile) allelefile="$2"; shift 2 ;;
        --reference) reference="$2"; shift 2 ;;
        -t|--threads) threads="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Error: Unknown option $1"; usage ;;
    esac
done

# Validate required arguments
if [ -z "$inputbam" ] || [ -z "$resultdir" ] || [ -z "$ref" ] || [ -z "$gtf" ] || [ -z "$kmer" ] || [ -z "$allelefile" ] || [ -z "$reference" ]; then
    echo "Error: All required arguments must be provided" >&2
    usage
fi

# validate --reference is actually hg19 or hg38
if [ "$reference" != "hg19" ] && [ "$reference" != "hg38" ]; then
    echo "Error: --reference must be exactly 'hg19' or 'hg38' (got: '$reference')" >&2
    exit 1
fi

# validate threads is a positive integer
if ! [[ "$threads" =~ ^[0-9]+$ ]] || [ "$threads" -lt 1 ]; then
    echo "Error: --threads must be a positive integer (got: '$threads')" >&2
    exit 1
fi

# Validate input files exist and are non-empty
for file in "$inputbam" "$ref" "$gtf" "$kmer" "$allelefile"; do
    if [ ! -f "$file" ]; then
        echo "Error: File not found: $file" >&2
        exit 1
    elif [ ! -s "$file" ]; then
        echo "Error: File is empty: $file" >&2
        exit 1
    fi
done

mkdir -p "$resultdir" || { echo "Error: could not create resultdir: $resultdir" >&2; exit 1; }

# check inputbam has an index before we do region-based samtools view later
if [ ! -f "${inputbam}.bai" ] && [ ! -f "${inputbam}.csi" ]; then
    echo "Warning: no .bai/.csi index found for $inputbam — indexing now" >&2

    if ! samtools index -@ "$threads" "$inputbam" 2>/dev/null; then
        rm -f "${inputbam}.bai" "${inputbam}.csi"   # clear any partial index
        echo "Warning: indexing failed (BAM not truly coordinate-sorted, header SO tag may be wrong) — sorting now" >&2

        sortedbam="${inputbam%.bam}.sorted.bam"
        samtools sort -@ "$threads" -T "${sortedbam%.bam}.tmp" -o "$sortedbam" "$inputbam" || {
            echo "Error: samtools sort failed for $inputbam" >&2; exit 1; }
        inputbam="$sortedbam"

        samtools index -@ "$threads" "$inputbam" || {
            echo "Error: samtools index failed for $inputbam even after sorting" >&2; exit 1; }
    fi
fi

echo "[1/9] Deduplicating allele list: $allelefile"
sort "$allelefile" | uniq > "${allelefile/.txt/.rnaseq.txt}"
sed -i '/^$/d' "${allelefile/.txt/.rnaseq.txt}"
allelefile="${allelefile/.txt/.rnaseq.txt}"

if [ ! -s "$allelefile" ]; then
    echo "Error: allele list is empty after dedup — nothing to align to" >&2
    exit 1
fi
n_alleles=$(wc -l < "$allelefile")
echo "    -> $n_alleles unique allele(s) to extract"

echo "[2/9] Extracting allele sequences into sample.fasta"
: > "$resultdir/sample.fasta"
fail_count=0
while IFS= read -r ALLELE; do
    if ! samtools faidx "$ref" "$ALLELE" >> "$resultdir/sample.fasta" 2>> "$resultdir/faidx.err"; then
        echo "Warning: could not extract allele '$ALLELE' from $ref (see $resultdir/faidx.err)" >&2
        fail_count=$((fail_count + 1))
    fi
done < "$allelefile"
if [ "$fail_count" -gt 0 ]; then
    echo "Warning: $fail_count/$n_alleles allele(s) failed extraction — check that allele names in $allelefile exactly match sequence names in $ref" >&2
fi
if [ ! -s "$resultdir/sample.fasta" ]; then
    echo "Error: sample.fasta is empty — no alleles were successfully extracted from $ref" >&2
    exit 1
fi

echo "[3/9] Building sample.gtf (filtering $gtf to requested alleles)"
: > "$resultdir/sample.gtf"
if ! grep -E -f "$allelefile" "$gtf" > "$resultdir/sample.gtf"; then
    echo "Error: no lines in $gtf matched any allele in $allelefile — check that allele naming conventions match between the two files" >&2
    exit 1
fi
if [ ! -s "$resultdir/sample.gtf" ]; then
    echo "Error: sample.gtf is empty after filtering" >&2
    exit 1
fi

echo "[4/9] Building STAR genome index (alleles only)"
readlen=$(samtools view -F 0x900 "$inputbam" | head -1000 | awk '$10!="*"{print length($10)}' | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')
sjdb_overhang=$((readlen - 1))

STAR --runMode genomeGenerate \
  --runThreadN "$threads" \
  --genomeDir "$resultdir/star_index" \
  --genomeFastaFiles "$resultdir/sample.fasta" \
  --sjdbGTFfile "$resultdir/sample.gtf" \
  --sjdbOverhang "$sjdb_overhang" \
  --genomeSAindexNbases 6 \
  || { echo "Error: STAR genomeGenerate failed — see $resultdir/star_index/Log.out" >&2; exit 1; }

echo "[5/9] Extracting kmer-matching reads"
samtools view -@ "$threads" -H "$inputbam" > "$resultdir/rna.kmer.sam"
samtools view -@ "$threads" "$inputbam" | grep -F -f "$kmer" >> "$resultdir/rna.kmer.sam"
grep_status=$?
if [ "$grep_status" -gt 1 ]; then
    echo "Error: grep failed while extracting kmer reads (exit code $grep_status)" >&2
    exit 1
elif [ "$grep_status" -eq 1 ]; then
    echo "Warning: zero reads matched the kmer file $kmer — continuing, but check this is expected" >&2
fi
samtools view -@ "$threads" -bS -o "$resultdir/rna.kmer.bam" "$resultdir/rna.kmer.sam" \
  || { echo "Error: failed to convert rna.kmer.sam to BAM" >&2; exit 1; }

awk -F'\t' '/^@/ || !seen[$1 FS $2]++' "$resultdir/rna.kmer.sam" \
  | samtools view -b -o "$resultdir/rna.kmer.dedup.bam" -

mv "$resultdir/rna.kmer.dedup.bam" "$resultdir/rna.kmer.bam" 

gatk SamToFastq -I "$resultdir/rna.kmer.bam" \
  -F "$resultdir/rna.kmer.1.fastq" -F2 "$resultdir/rna.kmer.2.fastq" \
  --VALIDATION_STRINGENCY SILENT

echo "[6/9] Extracting chr6 HLA-region reads (reference build: $reference)"
samtools view -@ "$threads" -H "$inputbam" > "$resultdir/rna.chr6region.sam"
chr_prefix=$(samtools view -H $inputbam | awk '/^@SQ/ && /SN:chr6/{found=1} END{print (found ? "chr" : "")}')
if [ $reference == "hg38" ]; then
    regions=("${chr_prefix}6:29941260-29945884" "${chr_prefix}6:31353872-31357187" "${chr_prefix}6:31268749-31272105")
elif [ $reference == "hg19" ]; then
    regions=("${chr_prefix}6:29909037-29913661" "${chr_prefix}6:31321649-31324964" "${chr_prefix}6:31236526-31239869")
fi
region_reads_total=0
for region in "${regions[@]}"; do
    n=$(samtools view -@ "$threads" "$inputbam" "$region" | tee -a "$resultdir/rna.chr6region.sam" | wc -l)
    if [ "$n" -eq 0 ]; then
        echo "Warning: zero reads found for region $region — check contig naming (chr6 vs 6) matches your BAM header" >&2
    fi
    region_reads_total=$((region_reads_total + n))
done
if [ "$region_reads_total" -eq 0 ]; then
    echo "Error: zero reads extracted across all chr6 HLA regions — likely a contig-naming mismatch (e.g. BAM uses '6' not 'chr6')" >&2
    exit 1
fi

samtools view -@ "$threads" -bS -o "$resultdir/rna.chr6region.bam" "$resultdir/rna.chr6region.sam" \
  || { echo "Error: failed to convert rna.chr6region.sam to BAM" >&2; exit 1; }

awk -F'\t' '/^@/ || !seen[$1 FS $2]++' "$resultdir/rna.chr6region.sam" \
  | samtools view -b -o "$resultdir/rna.chr6region.dedup.bam" -

mv "$resultdir/rna.chr6region.dedup.bam" "$resultdir/rna.chr6region.bam"

gatk SamToFastq -I "$resultdir/rna.chr6region.bam" \
  -F "$resultdir/rna.chr6region.1.fastq" -F2 "$resultdir/rna.chr6region.2.fastq" \
  --VALIDATION_STRINGENCY SILENT


echo "[7/9] Merging and deduplicating R1/R2 fastqs"
cat "$resultdir/rna.kmer.1.fastq" "$resultdir/rna.chr6region.1.fastq" | awk '{
    if(NR%4==1) {
      id=$1
      header=$0
      getline seq
      getline plus
      getline qual
      if(!seen[id]++) {
        print header"\n"seq"\n"plus"\n"qual
      }
    }
  }' > "$resultdir/rna.R1.fastq"

cat "$resultdir/rna.kmer.2.fastq" "$resultdir/rna.chr6region.2.fastq" | awk '{
    if(NR%4==1) {
      id=$1
      header=$0
      getline seq
      getline plus
      getline qual
      if(!seen[id]++) {
        print header"\n"seq"\n"plus"\n"qual
      }
    }
  }' > "$resultdir/rna.R2.fastq"

if [ ! -s "$resultdir/rna.R1.fastq" ] || [ ! -s "$resultdir/rna.R2.fastq" ]; then
    echo "Error: merged R1/R2 fastqs are empty — no reads survived kmer + chr6-region extraction" >&2
    exit 1
fi
r1_reads=$(( $(wc -l < "$resultdir/rna.R1.fastq") / 4 ))
echo "    -> $r1_reads read pairs going into STAR alignment"

echo "[8/9] Aligning to allele-specific STAR index"
run_star() {
    local prefix=$1; shift
    local retries=${STAR_SORT_RETRIES:-1}
    local headroom=${STAR_SORT_HEADROOM_PCT:-10}
    local attempt=0 status suggested newram arg skip
    local -a cmd=("$@") next

    while :; do
        "${cmd[@]}"
        status=$?
        [ "$status" -eq 0 ] && return 0

        suggested=$(grep -o 'limitBAMsortRAM [0-9]\+' "${prefix}Log.out" 2>/dev/null | tail -1 | awk '{print $2}')
        if [ -z "$suggested" ]; then
            echo "Error: STAR alignment failed (exit $status) — see ${prefix}Log.out" >&2
            return "$status"
        fi
        if [ "$attempt" -ge "$retries" ]; then
            echo "Error: STAR BAM sorting still out of memory after $attempt retry(ies) (last suggestion: $suggested bytes) — see ${prefix}Log.out" >&2
            return "$status"
        fi

        attempt=$((attempt + 1))
        newram=$(( suggested + suggested * headroom / 100 ))
        echo "Warning: STAR ran out of memory for BAM sorting — retry $attempt with --limitBAMsortRAM $newram" >&2

        next=(); skip=0
        for arg in "${cmd[@]}"; do
            if [ "$skip" -eq 1 ]; then skip=0; continue; fi
            if [ "$arg" = "--limitBAMsortRAM" ]; then skip=1; continue; fi
            next+=("$arg")
        done
        cmd=("${next[@]}" --limitBAMsortRAM "$newram")

        rm -rf "${prefix}_STARtmp" "${prefix}_STARpass1" "${prefix}_STARgenome"
    done
}

run_star "$resultdir/" STAR \
  --genomeDir "$resultdir/star_index" \
  --readFilesIn "$resultdir/rna.R1.fastq" "$resultdir/rna.R2.fastq" \
  --runThreadN "$threads" \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM GeneCounts \
  --twopassMode Basic \
  --sjdbGTFfile "$resultdir/sample.gtf" \
  --sjdbOverhang "$sjdb_overhang" \
  --outFileNamePrefix "$resultdir/" \
  || exit 1
if [ ! -s "$resultdir/Aligned.sortedByCoord.out.bam" ]; then
    echo "Error: STAR did not produce Aligned.sortedByCoord.out.bam" >&2
    exit 1
fi


echo "[9/9] Indexing and MAPQ-correcting multi-mappers"
samtools index -@ "$threads" "$resultdir/Aligned.sortedByCoord.out.bam"
samtools view -h "$resultdir/Aligned.sortedByCoord.out.bam" > "$resultdir/Aligned.sortedByCoord.out.sam"

awk 'BEGIN {OFS="\t"}
     /^@/ {print; next}
     {
       nh=0
       for(i=12; i<=NF; i++) {
         if($i ~ /^NH:i:/) {
           split($i, a, ":")
           nh=a[3]
           break
         }
       }
       if(nh > 1) $5=0  # Set MAPQ to 0 for multi-mappers
       print
     }' "$resultdir/Aligned.sortedByCoord.out.sam" | \
  samtools view -bS > "$resultdir/Aligned.sortedByCoord.out_mapqcorrected.bam"
samtools index -@ "$threads" "$resultdir/Aligned.sortedByCoord.out_mapqcorrected.bam"

samtools view -@ $threads -h "$resultdir/Aligned.sortedByCoord.out_mapqcorrected.bam" > "$resultdir/Aligned.sortedByCoord.out_mapqcorrected.sam"

echo "Creating separate bam files for MAPQzero and MAPQnonzero reads ..."
samtools view -@ $threads -H "$resultdir/Aligned.sortedByCoord.out_mapqcorrected.sam" > "$resultdir/rna.tumor.MAPQzero.sam"
grep -v ^@ "$resultdir/Aligned.sortedByCoord.out_mapqcorrected.sam" | awk '$5==0{print $0}' >> "$resultdir/rna.tumor.MAPQzero.sam"
samtools view -@ $threads -bS -o "$resultdir/rna.tumor.MAPQzero.bam" "$resultdir/rna.tumor.MAPQzero.sam"
samtools index -@ $threads "$resultdir/rna.tumor.MAPQzero.bam"

samtools view -@ $threads -H "$resultdir/Aligned.sortedByCoord.out_mapqcorrected.sam" > "$resultdir/rna.tumor.MAPQnonzero.sam"
grep -v ^@ "$resultdir/Aligned.sortedByCoord.out_mapqcorrected.sam" | awk '$5>0{print $0}' >> "$resultdir/rna.tumor.MAPQnonzero.sam"
samtools view -@ $threads -bS -o "$resultdir/rna.tumor.MAPQnonzero.bam" "$resultdir/rna.tumor.MAPQnonzero.sam"
samtools index -@ $threads "$resultdir/rna.tumor.MAPQnonzero.bam"


echo "Done. Final BAMs: $resultdir/Aligned.sortedByCoord.out_mapqcorrected.bam, $resultdir/rna.tumor.MAPQzero.bam, and $resultdir/rna.tumor.MAPQnonzero.bam"
