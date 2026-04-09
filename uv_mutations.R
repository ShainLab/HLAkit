#!/usr/bin/env Rscript
#!/usr/local/bin Rscript

library(optparse)
library(Rsamtools)

option_list <- list( 
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
mml <- opt$mml
fastafile <- opt$fastafile


if (is.null(mml)) stop("ERROR: -m/--mml is required but not provided.")
if (is.null(fastafile)) stop("ERROR: -f/--fastafile is required but not provided.")

check_file <- function(f) {
    if (!file.exists(f)) stop(paste("File not found:", f))
    if (file.info(f)$size == 0) stop(paste("File is empty:", f))
}
check_file(mml)
check_file(fastafile)

# Check FASTA index exists
fai <- paste0(fastafile, ".fai")
if (!file.exists(fai)) stop(paste("ERROR: FASTA index (.fai) not found:", fai, "— run 'samtools faidx' on the FASTA file before running this script."))

cat("Input files check complete.\n")

outfile <- sub('.txt', '.UV.txt', mml)
seqcontext_file <- sub('.txt', '.UV.seqcontext.txt', mml)

mml <- tryCatch(read.delim(mml, sep = "\t", header = T), error=function(e) NULL)
if (is.null(mml)) stop("ERROR: Failed to read MML file. Check that it is a valid tab-delimited file with a header.")
if (nrow(mml) == 0) stop("ERROR: MML file has no data rows.")

required_cols <- c("CHROM", "POS", "REF", "ALT")
missing_cols <- setdiff(required_cols, colnames(mml))
if (length(missing_cols) > 0) stop(paste("ERROR: MML file is missing required column(s):", paste(missing_cols, collapse = ", ")))

if (any(is.na(mml$POS))) {
  n_na <- sum(is.na(mml$POS))
  warning(paste0("WARNING: ", n_na, " row(s) have NA in the POS column. These will produce NA GRanges coordinates and may cause scanFa to fail or return empty sequence."))
}

gpos <- as.data.frame(cbind("chr" = paste0(mml[ , "CHROM"]), "pos1" = as.numeric(mml[ , "POS"]) - 1, "pos2" = as.numeric(mml[ , "POS"]) + 1))
gpos <- GenomicRanges::makeGRangesFromDataFrame(gpos, start.field = 'pos1', end.field = 'pos2')

missing_base <- tryCatch(
  as.data.frame(scanFa(fastafile, gpos)),
  error = function(e) stop(paste("ERROR: scanFa failed. Check that CHROM values in MML match sequence names in the FASTA file, and that the .fai index is up to date. Details:", conditionMessage(e)))
)

if (nrow(missing_base) == 0) stop("ERROR: scanFa returned no sequences. Check that CHROM values in the MML match sequence names in the FASTA.")
if (nrow(missing_base) != nrow(mml)) stop(paste0("ERROR: scanFa returned ", nrow(missing_base), " sequences but MML has ", nrow(mml), " rows. These must match."))

write.table(cbind(mml, missing_base), file=seqcontext_file, sep = '\t', quote=F, row.names =F)

missing_base <- as.data.frame(do.call('rbind', lapply(1:nrow(missing_base), function(x) as.vector(strsplit(missing_base[x, 1], split='', fixed=T))[[1]] )))

if (ncol(missing_base) < 3) stop(paste0("ERROR: Sequence context extraction produced only ", ncol(missing_base), " column(s); expected 3 (left base, ref base, right base). Check that POS values are not at the very start/end of the contig."))

index_c_t <- which(mml[,"REF"] == "C" & mml[,"ALT"] == "T")
index_g_a <- which(mml[,"REF"] == "G" & mml[,"ALT"] == "A")
index_cc_tt_gg_aa <- which( (mml[,"REF"] == "CC" & mml[,"ALT"] == "TT") | (mml[,"REF"] == "GG" & mml[,"ALT"] == "AA"))

message(paste0("NOTE: UV candidate counts — C>T: ", length(index_c_t),
               ", G>A: ", length(index_g_a),
               ", CC>TT/GG>AA: ", length(index_cc_tt_gg_aa)))

mml$UV <- NA
if(length(index_c_t) > 0){
  for(i in 1:length(index_c_t)){
    if( (missing_base[index_c_t[i],1] == "C") || (missing_base[index_c_t[i],1] == "T")) mml[index_c_t[i],'UV'] <- "UV"
  }
} 
if(length(index_g_a) > 0){
  for(i in 1:length(index_g_a)){
    if( (missing_base[index_g_a[i],3] == "G") || (missing_base[index_g_a[i],3] == "A")) mml[index_g_a[i],'UV'] <- "UV"
  }
}
if(length(index_cc_tt_gg_aa) > 0) mml[index_cc_tt_gg_aa, 'UV'] <- 'UV'
mml$UV[is.na(mml$UV)] <- ""

n_uv <- sum(mml$UV == "UV")
message(paste0("NOTE: ", n_uv, " / ", nrow(mml), " mutation(s) annotated as UV."))
if (n_uv == 0) warning("WARNING: No mutations were annotated as UV.")

cat("UV column added successfully. Writing output to ", outfile, "\n")
write.table(mml, sep = "\t", row.names = F, quote = F, file = outfile)

                                                      
