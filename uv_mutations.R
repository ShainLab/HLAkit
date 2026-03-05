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

check_file <- function(f) {
    if (!file.exists(f)) stop(paste("File not found:", f))
    if (file.info(f)$size == 0) stop(paste("File is empty:", f))
}

check_file(mml)
check_file(fastafile)


cat("Input files check complete.\n")


outfile <- sub('.txt', '.UV.txt', mml)
seqcontext_file <- sub('.txt', '.UV.seqcontext.txt', mml)
mml <- tryCatch(read.delim(mml, sep = "\t", header = T), error=function(e) NULL)

gpos <- as.data.frame(cbind("chr" = paste0(mml[ , "CHROM"]), "pos1" = as.numeric(mml[ , "POS"]) - 1, "pos2" = as.numeric(mml[ , "POS"]) + 1))
gpos <- GenomicRanges::makeGRangesFromDataFrame(gpos, start.field = 'pos1', end.field = 'pos2')

missing_base <- as.data.frame(scanFa(fastafile, gpos))
write.table(cbind(mml, missing_base), file=seqcontext_file, sep = '\t', quote=F, row.names =F)

missing_base <- as.data.frame(do.call('rbind', lapply(1:nrow(missing_base), function(x) as.vector(strsplit(missing_base[x, 1], split='', fixed=T))[[1]] )))


index_c_t <- which(mml[,"REF"] == "C" & mml[,"ALT"] == "T")
index_g_a <- which(mml[,"REF"] == "G" & mml[,"ALT"] == "A")
index_cc_tt_gg_aa <- which( (mml[,"REF"] == "CC" & mml[,"ALT"] == "TT") | (mml[,"REF"] == "GG" & mml[,"ALT"] == "AA"))

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

if(length(index_cc_tt_gg_aa) > 1) mml[index_cc_tt_gg_aa, 'UV'] <- 'UV'

mml$UV[is.na(mml$UV)] <- ""

cat("UV column added successfully. Writing output to ", outfile, "\n")
write.table(mml, sep = "\t", row.names = F, quote = F, file = outfile)

