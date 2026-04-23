# HLAkit
Software to detect both somatic point mutations and loss of heterozygosity affecting HLA class I genes

## Platform
Linux

## Contents
- R scripts for annotation, data analysis and generating plots
- Bash scripts for data analysis and pipeline automation

## Requirements
HLAkit needs the following bash tools and R packages to be installed:

### BASH tools:
The following tools should be installed and in PATH:
1. samtools (REQUIRED):
To install samtools, follow the instructions here: https://www.htslib.org/download . Also install the samtools dependencies listed on the link (bcftools and htslib)
2. gatk4 (REQUIRED):
To install GATK4, follow the instructions here: https://github.com/broadinstitute/gatk
4. novoalign (REQUIRED):
Download the novoalign binary from here: https://www.novocraft.com/support/download
6. R4.4+ (REQUIRED)
7. GNU parallel (OPTIONAL):
To install GNU parallel, follow the instructions here: https://www.gnu.org/software/parallel

### R packages:
1. dplyr (REQUIRED)
2. Rsamtools (REQUIRED)
3. stringr (REQUIRED)
4. readr (REQUIRED)
5. IRanges (REQUIRED)
6. tidyverse (REQUIRED)
7. Biostrings (REQUIRED)
8. GenomicRanges (REQUIRED)
9. optparse (REQUIRED)
10. tidyr (REQUIRED)
11. ggplot2 (REQUIRED)

## Installation:
Clone the github repo:
```
git clone https://github.com/ShainLab/HLAkit
```
Make all the scripts executable:
```
cd HLAkit
chmod +x hlakit *sh *R binaries/*
```
Uncompress the hla fasta file:
```
hlakit=$PWD
tar -xzvf $hlakit/resources/hla.fasta.tar.gz -C $hlakit/resources
```
To create fasta index for novoalign, download and install novoindex from here: https://www.novocraft.com/support/download
```
novoindex $hlakit/resources/hla.nix $hlakit/resources/hla.fasta
```

## Instructions for creating allelelist
Convert HLA allele names from Standard IMGT format (`A*01:01:01:01`) to flattened format (`hla_a_01_01_01_01`). Also remove the expression status suffix from allele name, for example `A*02:01:01:01N` -> `hla_a_02_01_01_01`.
```
hla_to_flat <- function(x) {
  x |>
    sub("[A-Za-z]$", "", x = _) |> 
    tolower() |>
    sub("\\*", "_", x = _) |>
    gsub(":", "_", x = _, fixed = TRUE) |>
    paste0("hla_", x = _)
}

hla_to_flat("A*02:01:01:01N")
# → "hla_a_02_01_01_01"
```

## Calculating coverage of exome-wide Normal and Tumor bams
Use gatk CollectHsMetrics to calculate exome-wide coverage and use MEAN_BAIT_COVERAGE as input for WESnormalcoverage and WEStumorcoverage.

## Grooming WES Normal and Tumor bam files
Before running HLAkit, make sure to deduplicate the WES bam files using `gatk MarkDuplicates` to reduce the run time and memory usage.

## Run:
Use 4 or 8 threads for multithreading.

hlakit=/path/to/hlakit/directory

Only somatic mutation calling:
```
$hlakit/hlakit --resultdir /path/to/resultdir --normalbam /path/to/normalbam --tumorbam /path/to/tumorbam --allelefile /path/to/allelelist --reference hg19 --hlakit $hlakit --tumorpurity tumorpurity_in_fraction --threads 8
```
Somatic mutation calling and Loss of Heterozygosity:
```
$hlakit/hlakit --resultdir /path/to/resultdir --normalbam /path/to/normalbam --tumorbam /path/to/tumorbam --allelefile /path/to/allelelist --reference hg19 --hlakit $hlakit --threads 8 --loh yes --tumorpurity tumorpurity_in_fraction --WESnormalcoverage WESnormalcoverage --WEStumorcoverage WEStumorcoverage
```

## Output files:
1. result/bamfiles/SNP
   Bam files used to find SNPs, Indels, and coverage.
   normal.SNP.bam: allele specific normal bam file with MAPQ=70 for all reads. This bam file is used to run Mutect2 and to count Normal Ref and Alt counts.
   tumor.SNP.bam: allele specific tumor bam file with MAPQ=70 for all reads. This bam file is used to run Mutect2 and to count Tumor Alt counts.
   normal.SNP.MAPQcorrected.bam: allele specific normal bam file with MAPQ=0 for multimapping reads and MAPQ=70 reads mapping uniquely to the allele. This bam file is split into MAPQ=0 and MAPQ=70 reads to find Normal allele-specific coverage.
   tumor.SNP.MAPQcorrected.bam: allele specific tumor bam file with MAPQ=0 for multimapping reads and MAPQ=70 reads mapping uniquely to the allele. This bam file is used to find Tumor Ref counts and tumor allele-specific coverage.
2. result/bamfiles/DNP
   Bam files used to find DNPs.
   normalDSNP.bam: allele specific normal bam file with MAPQ=70 for all reads.
   tumor.DNP.bam: allele specific tumor bam file with MAPQ=70 for all reads.
   normal.DNP.MAPQcorrected.bam: allele specific normal bam file with MAPQ=0 for multimapping reads and MAPQ=70 reads mapping uniquely to the allele.
   tumor.DNP.MAPQcorrected.bam: allele specific tumor bam file with MAPQ=0 for multimapping reads and MAPQ=70 reads mapping uniquely to the allele.
3. result/bamfiles/allelelist_updated.txt
   Nearest alleles with a genomic sequence if any input allele did not have a genomic sequence in hla.fasta
4. result/somatic_mutations/SNP
   Mutect2 output VCF files for each allele.
5. result/somatic_mutations/DNP
   Mutect2 output VCF files for each allele. HLAkit discards any SNP or Indel from these VCF files and only considers DNPs.
6. result/somatic_mutations/hla_somatic_mutations.txt
   Compiled mutations in all alleles of the sample with Mutect2 Filter, Normal/Tumor Ref and Alt counts, Tumor Mutant Allele Frequency (MAF), Gene Feature in which mutation is detected, along with annotations for - variant_type, variant_classification, UV-induced mutation, amino_acid_change, CDS_position, amino_acid_position, Zygosity of the gene, zygosity-corrected Tumor Mutant Allele Frequency, and Artifacts assigned by HLAkit based on Normal Alt count, Tumor Alt count, and Tumor Mutant Allele Frequency.
7. result/somatic_mutations/hla_somatic_mutations_filt.txt
   Filtered somatic mutations derived from somatic_mutations/hla_somatic_mutations.txt. Excludes mutations lacking a PASS filter flag or flagged as artifacts by HLAkit.
8. result/coverage/coverage.txt
   Allele-specific coverage in Normal and Tumor bam files.
9. result/loh/loh.txt
   Loss of heterozygosity result of HLAkit.
   aib_or_lowcov=NA  - HLAkit could not make a call due to low coverage
   aib_or_lowcov=0   - No LOH
   aib_or_lowcov=0.5 - Allelic imbalance (unequal normalized tumor coverage across alleles) or low coverage (normalized tumor coverage lower than normal) detected
   aib_or_lowcov=1   - Allelic imbalance (unequal normalized tumor coverage across alleles) and low coverage (normalized tumor coverage lower than normal) detected
   
   

## Test:

Only somatic mutation calling:
```
$hlakit/hlakit --resultdir $hlakit/test --normalbam $hlakit/test/normalWES.bam --tumorbam $hlakit/test/tumorWES.bam --allelefile $hlakit/test/allelelist.txt --reference hg19 --hlakit $hlakit --tumorpurity 0.513 --threads 8
```
Somatic mutation calling and Loss of Heterozygosity:
```
$hlakit/hlakit --resultdir $hlakit/test --normalbam $hlakit/test/normalWES.bam --tumorbam $hlakit/test/tumorWES.bam --allelefile $hlakit/test/allelelist.txt --reference hg19 --hlakit $hlakit --threads 8 --loh yes --tumorpurity 0.513 --WESnormalcoverage 166.054302 --WEStumorcoverage 289.227264
```

## Visualize the bam file in IGV:

To visualize bam files generated by HLAkit:
1. Update the absolute paths to hla.fasta, hla.fasta.fai, and hla_igv.bed in hlakit.igv.json file.
2. Open IGV -> Genomes -> Load Genome from File -> hlakit.igv.json.
3. Load the required bam file(s).
