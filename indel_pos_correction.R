#!/usr/bin/env Rscript
#!/usr/local/bin Rscript

# indel pos correction
library(optparse)
library(Rsamtools)

option_list <- list(
    make_option(c("-o", "--resultdir"), 
        type = "character",
        help = "result directory"
        ),
    make_option(c("-m", "--mml"), 
        type = "character",
        help = "somatic mutations file"
        ),
    make_option(c("-f", "--fastafile"), 
        type = "character",
        help = "fasta file"
        )
)

# read the arguments
opt <- parse_args(OptionParser(option_list = option_list))
# allot files to variables
resultdir <- opt$resultdir
mmlfile <- opt$mml
fastafile <- opt$fastafile


if (is.null(mmlfile)) stop("ERROR: -m/--mml is required but not provided.")
if (is.null(resultdir)) stop("ERROR: -o/--resultdir is required but not provided.")
if (is.null(fastafile)) stop("ERROR: -f/--fastafile is required but not provided.")

check_file <- function(f) {
    if (!file.exists(f)) stop(paste("File not found:", f))
    if (file.info(f)$size == 0) stop(paste("File is empty:", f))
}
check_file(mmlfile)
check_file(fastafile)

# Check FASTA index exists
fai <- paste0(fastafile, ".fai")
if (!file.exists(fai)) stop(paste("ERROR: FASTA index (.fai) not found:", fai, "— run 'samtools faidx' on the FASTA file before running this script."))

cat("Input files check complete.\n")


mml <- tryCatch(read.delim(mmlfile, sep = "\t", header = T, colClasses = c(REF = "character", ALT = "character")), error=function(e) NULL)
if (is.null(mml)) stop("ERROR: Failed to read MML file. Check that it is a valid tab-delimited file with a header.")
if (nrow(mml) == 0) stop("ERROR: MML file has no data rows.")

required_cols <- c("CHROM", "POS", "REF", "ALT")
missing_cols <- setdiff(required_cols, colnames(mml))
if (length(missing_cols) > 0) stop(paste("ERROR: MML file is missing required column(s):", paste(missing_cols, collapse = ", ")))

if (any(is.na(mml$POS))) {
  n_na <- sum(is.na(mml$POS))
  warning(paste0("WARNING: ", n_na, " row(s) have NA in the POS column. These will produce NA GRanges coordinates and may cause scanFa to fail or return empty sequence."))
}

# function to update indel position
get_true_indel <- function(bam_path, chrom, pos, window=10) {
  gr <- GRanges(chrom, IRanges(pos - window, pos + window))
  
  p <- pileup(bam_path,
              scanBamParam = ScanBamParam(which = gr),
              pileupParam  = PileupParam(distinguish_strands = FALSE,
                                         include_insertions  = TRUE))
  
  # insertions show up as "+" prefixed bases
  ins <- p[grepl("^\\+", p$nucleotide), ]
  del <- p[p$nucleotide == "-", ]
  
  indels <- rbind(ins, del)
  if (nrow(indels) == 0) { message("No indels found near ", chrom, ":", pos); return(NULL) }
  
  # position with most supporting reads
  best <- indels[which.max(indels$count), ]
  message("Reported pos: ", pos, " | True pos: ", best$pos, " | ", best$nucleotide, " (", best$count, " reads)")
  return(best)
}

# find correct ref base based on updated position
get_ref_base <- function(fasta, chrom, pos) {
    gr <- GRanges(seqnames = chrom, ranges = IRanges(start = pos, end = pos))
    return(as.character(scanFa(fasta, gr)))
}


# update position and first base of ref and alt in mml
if(nrow(mml) == 0){
  stop("No mutations found.")
} else{
  for(i in 1:nrow(mml)){
    allele <- mml$CHROM[i]
    ref <- mml$REF[i]
    alt <- mml$ALT[i]
    pos <- mml$POS[i]
    if(nchar(ref) != nchar(alt)){
      bamfile <- paste(resultdir, "/", allele, ".tumor.SNP.bam", sep = "", collapse="")
      if (!file.exists(bamfile)) stop(paste("BAM file not found:", bamfile))
      if (!file.exists(paste0(bamfile, ".bai"))) stop(paste("BAM index not found:", bamfile))
      indel_pos <- get_true_indel(bam_path = bamfile, chrom = allele, pos = pos)
      if (!is.null(indel_pos)) mml$POS[i] <- indel_pos$pos
      ref_base <- unname(get_ref_base(fastafile, allele, mml$POS[i]))
      substr(mml$REF[i], 1, 1) <- ref_base
      substr(mml$ALT[i], 1, 1) <- ref_base
    }
  }
}

write.table(mml, mmlfile, sep = "\t", quote = F, row.names = F)
