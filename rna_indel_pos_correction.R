#!/usr/bin/env Rscript

# Correct indel pos

suppressPackageStartupMessages({
    library(optparse)
    library(Rsamtools)
})

option_list <- list(
    make_option(c("-m", "--mutations"), type = "character", help = "somatic mutations file (not modified)"),
    make_option(c("-b", "--bam"),       type = "character", help = "RNA BAM to search for indels"),
    make_option(c("-o", "--outfile"),   type = "character", help = "output: mutation list with RNA-corrected indel POS"),
    make_option(c("-w", "--window"),    type = "integer", default = 10L, help = "search window +/- bp [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (a in c("mutations", "bam", "outfile")) {
    if (is.null(opt[[a]])) stop("ERROR: --", a, " is required but not provided.")
}
if (!file.exists(opt$mutations)) stop("ERROR: mutations file not found: ", opt$mutations)
if (!file.exists(opt$bam)) stop("ERROR: BAM not found: ", opt$bam)
if (!file.exists(paste0(opt$bam, ".bai"))) stop("ERROR: BAM index not found: ", opt$bam, ".bai")

mut <- read.delim(opt$mutations, sep = "\t", header = TRUE, quote = "",
                  colClasses = "character", check.names = FALSE)

missing_cols <- setdiff(c("CHROM", "POS", "REF", "ALT"), colnames(mut))
if (length(missing_cols) > 0) stop("ERROR: mutations file missing column(s): ", paste(missing_cols, collapse = ", "))

find_rna_indel_pos <- function(bam, chrom, pos, type, window) {
    gr <- GRanges(chrom, IRanges(max(1L, pos - window), pos + window))
    p <- pileup(bam,
                scanBamParam = ScanBamParam(which = gr),
                pileupParam  = PileupParam(max_depth = 1e6,
                                           min_base_quality = 0,
                                           min_mapq = 0,
                                           distinguish_strands = FALSE,
                                           include_deletions = TRUE,
                                           include_insertions = TRUE))

    # only consider indels of the same type as the reported mutation
    hits <- if (type == "insertion") p[grepl("^\\+", p$nucleotide), ] else p[p$nucleotide == "-", ]
    if (nrow(hits) == 0) return(NA_integer_)

    hits <- hits[order(-hits$count, abs(hits$pos - pos)), ]
    hits$pos[1]
}

is_ins <- nchar(mut$ALT) > nchar(mut$REF) | mut$REF %in% c("-", ".")
is_del <- nchar(mut$ALT) < nchar(mut$REF) | mut$ALT %in% c("-", ".")
indels <- which(is_ins | is_del)

message("NOTE: ", length(indels), " indel(s) to check for RNA position.")

for (i in indels) {
    chrom <- mut$CHROM[i]
    pos   <- as.integer(mut$POS[i])
    type  <- if (is_ins[i]) "insertion" else "deletion"

    if (is.na(pos)) {
        warning("Row ", i, ": POS is not numeric ('", mut$POS[i], "') - left unchanged.")
        next
    }

    new_pos <- tryCatch(find_rna_indel_pos(opt$bam, chrom, pos, type, opt$window),
                        error = function(e) {
                            warning("Row ", i, " (", chrom, ":", pos, "): pileup failed - ", conditionMessage(e))
                            NA_integer_
                        })

    if (is.na(new_pos)) {
        message(chrom, ":", pos, " ", type, " - no RNA support within +/-", opt$window, " bp, keeping original POS")
    } else {
        if (new_pos != pos) message(chrom, ":", pos, " ", type, " - RNA position ", new_pos)
        mut$POS[i] <- as.character(new_pos)
    }
}

write.table(mut, opt$outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cat("Indel position correction done:", opt$outfile, "\n")
