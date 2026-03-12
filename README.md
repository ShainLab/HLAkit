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

## Installation:
Clone the github repo:
```
git clone https://github.com/ShainLab/HLAkit
```
Make hlakit executable:
```
cd ShainLab/HLAkit
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
$hlakit/hlakit --resultdir /path/to/resultdir --normalbam /path/to/normalbam --tumorbam /path/to/tumorbam --allelefile /path/to/allelelist --reference hg19 --hlakit $hlakit --threads 8
```
Somatic mutation calling and Loss of Heterozygosity:
```
$hlakit/hlakit --resultdir /path/to/resultdir --normalbam /path/to/normalbam --tumorbam /path/to/tumorbam --allelefile /path/to/allelelist --reference hg19 --hlakit $hlakit --threads 8 --loh yes --tumorpurity tumorpurity_in_fraction --WESnormalcoverage WESnormalcoverage --WEStumorcoverage WEStumorcoverage
```

## Test:

Only somatic mutation calling:
```
$hlakit/hlakit --resultdir $hlakit/test --normalbam $hlakit/test/normalWES.bam --tumorbam $hlakit/test/tumorWES.bam --allelefile $hlakit/test/allelelist.txt --reference hg19 --hlakit $hlakit --threads 8
```
Somatic mutation calling and Loss of Heterozygosity:
```
$hlakit/hlakit --resultdir $hlakit/test --normalbam $hlakit/test/normalWES.bam --tumorbam $hlakit/test/tumorWES.bam --allelefile $hlakit/test/allelelist.txt --reference hg19 --hlakit $hlakit --threads 8 --loh yes --tumorpurity 0.513 --WESnormalcoverage 166.054302 --WEStumorcoverage 289.227264
```

