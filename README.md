# HLAkit
This repo contains code to find somatic mutations and loss of heterozygosity in HLA Class I genes.

## Contents
- R scripts for annotation, data analysis and generating plots
- Bash scripts for data analysis and pipeline automation

## Requirements
HLAkit needs the following bash tools and R packages to be installed:

### BASH tools:
1. samtools (REQUIRED)
To install samtools, follow the instructions here: https://www.htslib.org/download/
Also install the samtools dependencies listed on the link (bcftools and htslib)
2. gatk4 (REQUIRED)
To install GATK4, follow the instructions here: https://github.com/broadinstitute/gatk
4. novoalign (REQUIRED)
Download the novoalign binary here: https://www.novocraft.com/support/download
6. R4.4+ (REQUIRED)
7. GNU parallel (OPTIONAL)
To install GNU parallel, follow the instructions here: https://www.gnu.org/software/parallel/

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

# Installation:
Clone the github repo:
```
git clone https://github.com/ShainLab/HLAkit
cd ShainLab/HLAkit
chmod +x hlakit
```

# Run:
hlakit --resultdir outputdirectory --normalbam normalWES.bam --tumorbam tumorWES.bam --allelefile allelefile.txt --reference hg19 --hlakit /path/to/hlakit --format ILM1.8 --threads 8
