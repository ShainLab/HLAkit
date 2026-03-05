# HLAkit
This repo contains code to find somatic mutations and loss of heterozygosity in HLA Class I genes.

## Contents
- R scripts for annotation, data analysis and generating plots
- Bash scripts for data analysis and pipeline automation

## Requirements
HLAkit needs the following bash tools and R packages to be installed:

### BASH tools:
### samtools (REQUIRED)
To install samtools, follow the instructions here: https://www.htslib.org/download/
Also install the samtools dependencies listed on the link (bcftools and htslib)
### gatk4 (REQUIRED)
To install GATK4, follow the instructions here: https://github.com/broadinstitute/gatk
### novoalign (REQUIRED)
### GNU parallel (OPTIONAL)
To install GNU parallel, follow the instructions here: https://www.gnu.org/software/parallel/

### R packages:
#### R4.4+ (REQUIRED)
#### Rsamtools (REQUIRED)
#### stringr (REQUIRED)
#### readr (REQUIRED)
#### IRanges (REQUIRED)
#### tidyverse (REQUIRED)
#### Biostrings (REQUIRED)
#### GenomicRanges (REQUIRED)
#### optparse (REQUIRED)
#### tidyr (REQUIRED)
#### dplyr (REQUIRED)

# Installation:
Clone the github repo:
```
git clone https://github.com/ShainLab/HLAkit
cd ShainLab/HLAkit
chmod +x hlakit
```

# Run:
hlakit --resultdir outputdirectory --normalbam normalWES.bam --tumorbam tumorWES.bam --allelefile allelefile.txt --reference hg19 --hlakit /path/to/hlakit --format ILM1.8 --threads 8
