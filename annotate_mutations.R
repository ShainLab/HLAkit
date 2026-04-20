#!/usr/bin/env Rscript
#!/opt/homebrew/bin Rscript
#!/opt/homebrew/share/man/man1 Rscript.1

# Annotate MML with protein change, variant classification, CDS pos, amino acid pos, and variant type
library(optparse)
library(stringr)
library(readr)
library(tidyr)
library(dplyr)
library(IRanges)
library(tidyverse)
library(Biostrings)
library(GenomicRanges)



option_list <- list( 
    make_option(c("-m", "--mml"), 
        type = "character",
        help = "MAF file"
        ),
    make_option(c("-g", "--gtf"), 
        type = "character",
        help = "HLA annotation GTF file"
        ),
    make_option(c("-f", "--fastafile"), 
        type = "character",
        help = "HLA fasta file"
    ),
    make_option(c("-b", "--bed"),
    	type = "character",
    	help = "HLA annotation bed file"
    	),
    make_option(c("-p", "--tumorpurity"),
    	type = "numeric",
    	help = "Tumor purity in fraction")
)

# read the arguments
opt <- parse_args(OptionParser(option_list = option_list))


# allot files to variables
mml <- opt$mml
gtf <- opt$gtf
fastafile <- opt$fastafile
bed <- opt$bed
tumorpurity <- opt$tumorpurity

# --- Argument validation ---
if (is.null(mml))       		stop("ERROR: -m/--mml is required but not provided.")
if (is.null(gtf))       		stop("ERROR: -g/--gtf is required but not provided.")
if (is.null(fastafile)) 		stop("ERROR: -f/--fastafile is required but not provided.")
if (is.null(bed))       		stop("ERROR: -b/--bed is required but not provided.")
if (is.null(tumorpurity))   stop("ERROR: -b/--tumorpurity is required but not provided.")

check_file <- function(f) {
    if (!file.exists(f)) stop(paste("File not found:", f))
    if (file.info(f)$size == 0) stop(paste("File is empty:", f))
}

check_file(mml)
check_file(gtf)
check_file(fastafile)
check_file(bed)

# read files
outfile <- sub(".txt", ".annotated.txt", mml)
mml <- read.delim(mml, header = T, sep = '\t', quote = "", stringsAsFactors = F)
gtf <- read.delim(gtf, header = F, sep = '\t', quote = "", stringsAsFactors = F)
colnames(gtf) <- c("CHROM", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
genome <- readDNAStringSet(fastafile)
names(genome) <- sub(" .*", "", names(genome))  # Clean FASTA headers
bed <- read.delim(bed, header=F)
colnames(bed)[1:12] <- c("chrom", "start", "end", "name", "score", "strand",
                         "thickStart", "thickEnd", "itemRgb", "blockCount",
                         "blockSizes", "blockStarts")

# --- Post-read structural checks ---
if (nrow(mml) == 0) stop(paste("ERROR: MML file has no data rows:", opt$mml))
if (nrow(gtf) == 0) stop(paste("ERROR: GTF file has no data rows:", opt$gtf))
if (nrow(bed) == 0) stop(paste("ERROR: BED file has no data rows:", opt$bed))
if (length(genome) == 0) stop(paste("ERROR: FASTA file contains no sequences:", opt$fastafile))

required_mml_cols <- c("CHROM", "POS", "REF", "ALT")
missing_mml_cols <- setdiff(required_mml_cols, colnames(mml))
if (length(missing_mml_cols) > 0) stop(paste("ERROR: MML file is missing required column(s):", paste(missing_mml_cols, collapse = ", ")))

required_bed_cols <- 12
if (ncol(bed) < required_bed_cols) stop(paste0("ERROR: BED file has only ", ncol(bed), " column(s); expected at least 12 (standard BED12 format)."))

# Extract gene_id from attributes
gtf <- gtf %>%
filter(feature == "exon") %>%
mutate(gene_id = str_extract(attribute, 'gene_id "(.*?)"')) %>%
mutate(gene_id = str_remove_all(gene_id, 'gene_id |"')) 

if (nrow(gtf) == 0) stop("ERROR: No exon features found in GTF file after filtering. Check that the GTF contains rows with feature == 'exon'.")
if (all(is.na(gtf$gene_id))) stop("ERROR: gene_id could not be extracted from any GTF attribute field. Check that attributes follow the format: gene_id \"<id>\".")

# --- Check CHROM values in mml have matching gene_ids in GTF ---
mml_chroms <- unique(mml$CHROM)
gtf_genes  <- unique(gtf$gene_id)
missing_chroms <- setdiff(mml_chroms, gtf_genes)
if (length(missing_chroms) > 0) {
  warning(paste("WARNING: The following CHROM value(s) in the MML file have no matching gene_id in the GTF and will return 'unknown' annotations:",
                paste(missing_chroms, collapse = ", ")))
}

# --- Check CHROM values in mml have matching sequences in FASTA ---
missing_fasta <- setdiff(mml_chroms, names(genome))
if (length(missing_fasta) > 0) {
  warning(paste("WARNING: The following CHROM value(s) in the MML file have no matching sequence in the FASTA file:",
                paste(missing_fasta, collapse = ", ")))
}


dnp <- function(exons, pos, transcript_start, transcript_end, alt, cds){
	cds_pos <- 0
	genomic_to_cds <- list()
	for (i in 1:nrow(exons)) {
		for (p in exons$start[i]:exons$end[i]) {
			cds_pos <- cds_pos + 1
			genomic_to_cds[[as.character(p)]] <- cds_pos
		}
	}
	  
	# Determine codon affected
	in_exon_pos1 <- any(pos >= exons$start & pos <= exons$end) 
	in_exon_pos2 <- any((pos+1) >= exons$start & (pos+1) <= exons$end)
	utr_pos1 <- pos < transcript_start || pos > transcript_end
	utr_pos2 <- (pos+1) < transcript_start || (pos+1) > transcript_end
	splice_site <- any(abs(pos - exons$start) <= 2 | abs(pos - exons$end) <= 2 |
		abs((pos+1) - exons$start) <= 2 | abs((pos+1) - exons$end) <= 2)

	# splice site in introns and utrs
	if(splice_site && (!in_exon_pos1 || utr_pos1) && (!in_exon_pos2 || utr_pos2) )
		return(list(NA, NA, "splice_site", NA))

	# splice site; one base is in exon and the other in intron or utr
	if (splice_site &&
    (((!in_exon_pos1 || utr_pos1) && (in_exon_pos2 && !utr_pos2)) ||
     ((in_exon_pos1 && !utr_pos1) && (!in_exon_pos2 || utr_pos2)))) {
		if(in_exon_pos1) pos <- pos else pos <- pos+1
		codon_index <- floor((genomic_to_cds[[as.character(pos)]] - 1) / 3) + 1
		codon_start <- (codon_index - 1) * 3 + 1

		if (codon_start + 2 > length(cds)) return(list(NA, NA, "unknown", cds_pos))

		ref_codon <- as.character(subseq(cds, start = codon_start, width = 3))
		mut_cds <- cds
		cds_pos <- genomic_to_cds[[as.character(pos)]]
		mut_cds[cds_pos] <- DNAString(alt)
		mut_codon <- as.character(subseq(mut_cds, start = codon_start, width = 3))

		# Translate
		ref_aa <- GENETIC_CODE[ref_codon]
		mut_aa <- GENETIC_CODE[mut_codon]
		return(list(ref_aa, mut_aa, "splice_site", cds_pos))
	}

	# Intronic check
	if (!splice_site && !in_exon_pos1 && !in_exon_pos2) {
		return(list(NA, NA, "", NA))
	}

	# UTR check
	if (!splice_site && utr_pos1 && utr_pos2) {
		return(list(NA, NA, "utr", cds_pos))
	}

	cds_pos1 <- genomic_to_cds[[as.character(pos)]]
	cds_pos2 <- genomic_to_cds[[as.character(pos+1)]]

	if (!(as.character(pos) %in% names(genomic_to_cds))) return(list(NA, NA, "unknown", cds_pos))

	codon1_index <- floor((genomic_to_cds[[as.character(pos)]] - 1) / 3) + 1
	codon1_start <- (codon1_index - 1) * 3 + 1
	codon2_index <- floor((genomic_to_cds[[as.character(pos+1)]] - 1) / 3) + 1
	codon2_start <- (codon2_index - 1) * 3 + 1

	if (codon1_start + 2 > length(cds)) return(list(NA, NA, "unknown", cds_pos))

	ref_codon1 <- as.character(subseq(cds, start = codon1_start, width = 3))
	ref_codon2 <- as.character(subseq(cds, start = codon2_start, width = 3))
	mut_cds <- cds
	mut_cds[cds_pos1:cds_pos2] <- DNAString(alt)
  
	# if both are in same codon
	if(codon1_index == codon2_index){
	    mut_codon <- as.character(subseq(mut_cds, start = codon1_start, width = 3))
	    # Translate
	    ref_aa <- GENETIC_CODE[ref_codon1]
	    mut_aa <- GENETIC_CODE[mut_codon]

	    # Splice site check (±2 bp of any exon boundary)
		if (splice_site && (in_exon_pos1 & in_exon_pos2 & !utr_pos1 & !utr_pos2)){
		 	return(list(ref_aa, mut_aa, "splice_site", cds_pos))	 
		}

	    
	    if (ref_aa == mut_aa) {
	      return(list(ref_aa, mut_aa, "synonymous", cds_pos))
	      } else if (mut_aa == "*") {
	        return(list(ref_aa, mut_aa, "nonsense", cds_pos))
	        } else {
	          return(list(ref_aa, mut_aa, "missense", cds_pos))
	    	}
	}
	else{
	# if both are in different codons
	mut_codon1 <- as.character(subseq(mut_cds, start = codon1_start, width = 3))
	mut_codon2 <- as.character(subseq(mut_cds, start = codon2_start, width = 3))

	# Translate
	ref_aa1 <- GENETIC_CODE[ref_codon1]
	ref_aa2 <- GENETIC_CODE[ref_codon2]
	mut_aa1 <- GENETIC_CODE[mut_codon1]
	mut_aa2 <- GENETIC_CODE[mut_codon2]
	ref_aa <- paste(ref_aa1, ref_aa2, sep = "")
	mut_aa <- paste(mut_aa1, mut_aa2, sep = "")

	# Splice site check (±2 bp of any exon boundary)
	if (splice_site){
	 	return(list(ref_aa, mut_aa, "splice_site", cds_pos))	  	
	}

	if (ref_aa1 == mut_aa1 && ref_aa2 == mut_aa2) {
	  return(list(ref_aa, mut_aa, "synonymous", cds_pos))
	  } else if (mut_aa1 == "*" || mut_aa2 == "*") {
	    return(list(ref_aa, mut_aa, "nonsense", cds_pos))
	    } else {
	      return(list(ref_aa, mut_aa, "missense", cds_pos))
	    }
	  }

}

snp <- function(exons, pos, transcript_start, transcript_end, alt, cds){
	cds_pos <- 0
	in_exon <- any(pos >= exons$start & pos <= exons$end)
	utr <- pos < transcript_start || pos > transcript_end
	splice_site <- any(abs(pos - exons$start) <= 2 | abs(pos - exons$end) <= 2)


	# Intronic check
	if (!splice_site && !in_exon) {
		return(list(NA, NA, "", NA))
	}

	# UTR check
	if (!splice_site && utr) {
		return(list(NA, NA, "utr", NA))
	}

	# Splice site in intron or UTR
	if (splice_site && (utr || !in_exon)) {
		return(list(NA, NA, "splice_site", NA))
	}


	# Determine codon affected
	cds_pos <- 0
	genomic_to_cds <- list()
	for (i in 1:nrow(exons)) {
		for (p in exons$start[i]:exons$end[i]) {
			cds_pos <- cds_pos + 1
			genomic_to_cds[[as.character(p)]] <- cds_pos
		}
	}

	if (!(as.character(pos) %in% names(genomic_to_cds))) return(list(NA, NA, "unknown", cds_pos))

	codon_index <- floor((genomic_to_cds[[as.character(pos)]] - 1) / 3) + 1
	codon_start <- (codon_index - 1) * 3 + 1

	if (codon_start + 2 > length(cds)) return(list(NA, NA, "unknown", cds_pos))

	ref_codon <- as.character(subseq(cds, start = codon_start, width = 3))
	mut_cds <- cds
	cds_pos <- genomic_to_cds[[as.character(pos)]]
	mut_cds[cds_pos] <- DNAString(alt)
	mut_codon <- as.character(subseq(mut_cds, start = codon_start, width = 3))

	# Translate
	ref_aa <- GENETIC_CODE[ref_codon]
	mut_aa <- GENETIC_CODE[mut_codon]

	# Splice site check (±2 bp of any exon boundary)
	if (splice_site && (in_exon || !utr)){
		return(list(ref_aa, mut_aa, "splice_site", cds_pos))
	}


	if (ref_aa == mut_aa) {
		return(list(ref_aa, mut_aa, "synonymous", cds_pos))
		} else if (mut_aa == "*") {
			return(list(ref_aa, mut_aa, "nonsense", cds_pos))
			} else {
				return(list(ref_aa, mut_aa, "missense", cds_pos))
			}
}


annotate_variant <- function(chrom, pos, ref, alt) {
	# handling indels
	if(nchar(ref) != nchar(alt) || ref == "-" || alt == "-" ){
		if(nchar(ref) > nchar(alt) && nchar(ref) %% 3 == 0) return(list(NA, NA, "in_frame_deletion", NA))
		else if(nchar(ref) > nchar(alt) && nchar(ref) %% 3 != 0) return(list(NA, NA, "frame_shift_deletion", NA))
		else if(nchar(ref) < nchar(alt) && nchar(ref) %% 3 == 0) return(list(NA, NA, "in_frame_insertion", NA))
		else if(nchar(ref) < nchar(alt) && nchar(ref) %% 3 != 0) return(list(NA, NA, "frame_shift_insertion", NA))
		else return(list(NA, NA, "manually_inspect", NA))
	}

	gene_match <- gtf[which(gtf$gene_id == chrom),]

	if (nrow(gene_match) == 0) {
		warning(paste0("WARNING: No GTF exon records found for CHROM '", chrom, "' at POS ", pos, ". Returning 'unknown'."))
		return(list(NA, NA, "unknown", NA))
	}

	# Proceed with arranging
	exons <- gene_match %>% arrange(start)
	seqnames <- unique(exons$CHROM)
	valid_seqnames <- seqnames[seqnames %in% names(genome)]
	if (length(valid_seqnames) == 0) {
		warning(paste0("WARNING: No FASTA sequence found for any seqname (", paste(seqnames, collapse = ", "),
		               ") associated with CHROM '", chrom, "'. Returning 'unknown'."))
		return(list(NA, NA, "unknown", NA))
	}
	seqname <- valid_seqnames[1]

	if (!(seqname %in% names(genome))) return(list(NA, NA, "unknown", NA))

	gene_seq <- genome[[seqname]]
	strand <- unique(exons$strand)

	if (length(strand) > 1) {
		warning(paste0("WARNING: Multiple strand values (", paste(strand, collapse = ", "),
		               ") found for CHROM '", chrom, "'. Using first value: '", strand[1], "'."))
		strand <- strand[1]
	}

	  # Determine transcript boundaries
	  transcript_start <- min(exons$start)
	  transcript_end <- max(exons$end)

	  # Build CDS
	  cds <- DNAString()
	  for (i in 1:nrow(exons)) {
	  	cds <- append(cds, subseq(gene_seq, start = exons$start[i], end = exons$end[i]))
	  }
	  if (strand == "-") cds <- reverseComplement(cds)

	  if (length(cds) == 0) {
	  	warning(paste0("WARNING: CDS built for CHROM '", chrom, "' is empty. Check exon coordinates in GTF vs FASTA sequence length."))
	  	return(list(NA, NA, "unknown", NA))
	  }

########### DNP ###########
	is_dnp <- nchar(ref) == 2 && nchar(alt) == 2
	if (is_dnp) {
		return(dnp(exons=exons, pos=pos, transcript_start=transcript_start, transcript_end=transcript_end, alt=alt, cds=cds))
	}

########### SNP ###########
	is_snp <- nchar(ref) == 1 && nchar(alt) == 1
	  if(is_snp){
	  	return(snp(exons=exons, pos=pos, transcript_start=transcript_start, transcript_end=transcript_end, alt=alt, cds=cds))
	}

	warning(paste0("WARNING: Variant at CHROM '", chrom, "' POS ", pos, " REF '", ref, "' ALT '", alt,
	               "' did not match any handled variant type (SNP/DNP/indel). Returning 'manually_inspect'."))
	return(list(NA, NA, "manually_inspect", NA))
 }

result <- pmap(
	mml[, c("CHROM", "POS", "REF", "ALT")],
	function(CHROM, POS, REF, ALT) {
		annotate_variant(CHROM, POS, REF, ALT)
	}
)

result <- as.data.frame(do.call('rbind', result))
result[] <- lapply(result, as.character)
result[result == "NA"] <- ""
mml$amino_acid_change <- as.vector(paste(result[,1], result[,2], sep = ">"))
mml$amino_acid_change[mml$amino_acid_change == ">"] <- ""
mml$variant_classification <- result[,3]
mml$CDS_pos <- result[,4]


## amino acid position annotation
# Expand exon information
bed_expanded <- bed %>%
  mutate(
    blockSizes = strsplit(blockSizes, ","),
    blockStarts = strsplit(blockStarts, ",")
  ) %>%
  unnest(cols = c(blockSizes, blockStarts)) %>%
  mutate(
    blockSizes = as.numeric(blockSizes),
    blockStarts = as.numeric(blockStarts),
    exon_start = start + blockStarts,
    exon_end = exon_start + blockSizes
  )

if (nrow(bed_expanded) == 0) stop("ERROR: BED file produced no rows after expanding blockSizes/blockStarts. Check BED12 format.")

# Filter to coding exons only
coding_exons <- bed_expanded %>%
  filter(exon_end > thickStart & exon_start < thickEnd) %>%
  mutate(
    cds_start = pmax(exon_start, thickStart),
    cds_end = pmin(exon_end, thickEnd)
  )

if (nrow(coding_exons) == 0) stop("ERROR: No coding exons found after filtering BED file by thickStart/thickEnd. Check that thickStart and thickEnd define a valid CDS in the BED file.")

# Create GRanges for coding exons
exon_gr <- GRanges(
  seqnames = coding_exons$chrom,
  ranges = IRanges(start = coding_exons$cds_start, end = coding_exons$cds_end),
  exon_len = coding_exons$cds_end - coding_exons$cds_start,
  allele = coding_exons$chrom
)

# Convert mutations to GRanges
mut_gr <- GRanges(
  seqnames = mml$CHROM,
  ranges = IRanges(start = mml$POS, width = 1),
  strand = "*"
)

# Find overlaps
hits <- findOverlaps(mut_gr, exon_gr)

if (length(hits) == 0) {
  warning("WARNING: No overlaps found between mutations (mut_gr) and coding exons (exon_gr). All aa_pos values will be NA. Check that CHROM values in the MML match chrom values in the BED file, and that POS values fall within coding exon coordinates.")
}

# Prepare result
mml$aa_pos <- NA

for (i in seq_along(hits)) {
  mut_idx <- queryHits(hits)[i]
  exon_idx <- subjectHits(hits)[i]

  mut_start <- start(mut_gr[mut_idx])
  exon_start <- start(exon_gr[exon_idx])

  # Calculate nucleotide offset from coding start
  offset <- 0
  allele_exons <- exon_gr[as.character(seqnames(exon_gr)) == as.character(seqnames(mut_gr[mut_idx]))]

  # Sort by start
  allele_exons <- allele_exons[order(start(allele_exons))]

  for (j in seq_along(allele_exons)) {
    exon <- allele_exons[j]
    if (start(mut_gr[mut_idx]) >= start(exon) && start(mut_gr[mut_idx]) <= end(exon)) {
      # Inside this exon
      offset <- offset + (start(mut_gr[mut_idx]) - start(exon)) + 1
      break
    } else {
      # Add full exon length
      offset <- offset + width(exon)
    }
  }

  # Translate nucleotide offset to amino acid position
  aa_pos <- floor((offset + 2) / 3)
  mml$aa_pos[mut_idx] <- aa_pos - 1
}


mml <- as.data.frame(mml)
mml[mml == "NA"] <- ""


# add variant_type
mml$variant_type <- NA
mml$variant_type <- ifelse(nchar(mml$REF) == nchar(mml$ALT) & nchar(mml$ALT) == 1, "SNP",
                    ifelse(nchar(mml$REF) == nchar(mml$ALT) & nchar(mml$ALT) == 2, "DNP",
                    ifelse(nchar(mml$REF) > nchar(mml$ALT), "DEL",
                    ifelse(nchar(mml$REF) < nchar(mml$ALT), "INS", NA))))


## update variant_classification of indels in UTR and Introns
index <- which((mml$variant_type == "INS" | mml$variant_type == "DEL") & (mml$Feature == "utr" | grepl("^intron", mml$Feature)))
if(length(index) > 0) mml$variant_classification[index] <- NA

mml$Artifacts <- NA
mml$Artifacts <- ifelse(mml$Normal_Mut > 2, "artifact_high_normalMut",
    ifelse(mml$Tumor_Mut < 4 & mml$variant_type == "SNP", "artifact_low_tumorMut",
    	ifelse(mml$Tumor_Mut < 3 & mml$variant_type != "SNP", "artifact_low_tumorMut",
        ifelse(mml$Tumor_Mut + mml$Tumor_Ref == 0, "artifact_no_coverage",
        	ifelse(mml$Tumor_MAF < 0.4*tumorpurity, "artifact_low_tumorMAF", ""
        		)))))

# save output
write.table(mml, file = outfile, sep = "\t", quote = F, row.names = F, col.names = T)


cat("Annotation Done!\n")

