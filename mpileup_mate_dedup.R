#!/usr/bin/env Rscript
#!/opt/homebrew/bin Rscript
#!/opt/homebrew/share/man/man1 Rscript.1

#deduplicate overlapping mates
library(optparse)
library(Rsamtools)

option_list <- list( 
    make_option(c("-n", "--normal_mpileupout"), 
        type = "character",
        help = "Normal mpileup output"
        ),
    make_option(c("-m", "--tumor_mpileupout_mapq0"), 
        type = "character",
        help = "Tumor mpileup output for MAPQ0 bam"
    ),
    make_option(c("-u", "--tumor_mpileupout_mapqnonzero"),
        type = "character",
        help = "Tumor mpileup output for MAPQ nonzero bam")
)

# read the arguments
opt <- parse_args(OptionParser(option_list = option_list))


# allot files to variables
normal_mpileupout_file <- opt$normal_mpileupout
tumor_mpileupout_mapq0_file <- opt$tumor_mpileupout_mapq0
tumor_mpileupout_mapqnonzero_file <- opt$tumor_mpileupout_mapqnonzero

# --- Argument validation ---
if (is.null(normal_mpileupout_file)) stop("ERROR: -n/--normal_mpileupout_file is required but not provided.")
if (is.null(tumor_mpileupout_mapq0_file)) stop("ERROR: -m/--tumor_mpileupout_mapq0_file is required but not provided.")
if (is.null(tumor_mpileupout_mapqnonzero_file)) stop("ERROR: -u/--tumor_mpileupout_mapqnonzero_file is required but not provided.")

check_file <- function(f) {
    if (!file.exists(f)) stop(paste("File not found:", f))
    if (file.info(f)$size == 0) stop(paste("File is empty:", f))
}

check_file(normal_mpileupout_file)
check_file(tumor_mpileupout_mapq0_file)
check_file(tumor_mpileupout_mapqnonzero_file)

cat("Input files check complete.\n")

# read files
normal_mpileupout <- read.delim(normal_mpileupout_file, header = F, sep = '\t', quote = "", stringsAsFactors = F)
tumor_mpileupout_mapq0 <- read.delim(tumor_mpileupout_mapq0_file, header = F, sep = '\t', quote = "", stringsAsFactors = F)
tumor_mpileupout_mapqnonzero <- read.delim(tumor_mpileupout_mapqnonzero_file, header = F, sep = '\t', quote = "", stringsAsFactors = F)


if (nrow(normal_mpileupout) == 0) stop("ERROR: Normal mpileup output file has no data rows.")
if (nrow(tumor_mpileupout_mapq0) == 0) stop("ERROR: Tumor MAPQ0 mpileup output file has no data rows.")
if (nrow(tumor_mpileupout_mapqnonzero) == 0) stop("ERROR: Tumor MAPQ-nonzero mpileup output file has no data rows.")

#all mpileup files must match somatic mutations
# check mpileup files columns
for (nm in c("normal_mpileupout", "tumor_mpileupout_mapq0", "tumor_mpileupout_mapqnonzero")) {
  df <- get(nm)
  if (ncol(df) < 6) stop(paste0("ERROR: ", nm, " has only ", ncol(df), " column(s); expected at least 6 (chrom, pos, ref, depth, bases, qname). Check mpileup output format."))
}


parse_pileup_noqual <- function(bases, qnames) {
  qn <- strsplit(qnames, ",", fixed = TRUE)[[1]]
  obs <- vector("list", length(qn))
  i <- 1        # index into pileup string
  obs_i <- 1    # observation index
  while (i <= nchar(bases)) {
    ch <- substr(bases, i, i)
      # start-of-read marker ^ plus MAPQ char
      if (ch == "^") {
        i <- i + 2
        next
      }
        # end-of-read marker
        if (ch == "$") {
          i <- i + 1
          next
        }
        # deletion placeholder
        if (ch == "*") {
          obs[[obs_i]] <- list(
            allele = "*",
            qname  = qn[obs_i]
            )
          obs_i <- obs_i + 1
          i <- i + 1
          next
        }
        # ordinary base observation
        allele <- ch
        i <- i + 1
        # attached insertion/deletion
        if (i <= nchar(bases)) {
          ch2 <- substr(bases, i, i)
          if (ch2 %in% c("+", "-")) {
            sign <- ch2
            i <- i + 1
                # parse indel length
                lenstr <- ""
                while (
                  i <= nchar(bases) &&
                  grepl("[0-9]", substr(bases, i, i))
                  ) {
                  lenstr <- paste0(lenstr, substr(bases, i, i))
                  i <- i + 1
                }
                len <- as.integer(lenstr)
                indelseq <- substr(bases, i, i + len - 1)
                allele <- paste0(
                  allele,
                  sign,
                  len,
                  indelseq
                  )
                i <- i + len
              }
            }
            obs[[obs_i]] <- list(
              allele = allele,
              qname  = qn[obs_i]
              )
            obs_i <- obs_i + 1
          }
          obs
        }

dedup_overlap <- function(bases, qnames) {
  obs <- parse_pileup_noqual(bases, qnames)
  # keep first occurrence of each read name
  keep <- !duplicated(sapply(obs, `[[`, "qname"))
    obs2 <- obs[keep]
    data.frame(
    qname  = sapply(obs2, `[[`, "qname"),
    allele = sapply(obs2, `[[`, "allele"),
    stringsAsFactors = FALSE
  )
}


mpileupout_files <- c(normal_mpileupout_file, tumor_mpileupout_mapq0_file, tumor_mpileupout_mapqnonzero_file)
mpileupouts <- list(normal_mpileupout, tumor_mpileupout_mapq0, tumor_mpileupout_mapqnonzero)

for(i in 1:3){
  mpileupout <- mpileupouts[[i]]
  mpileupoutfile <- mpileupout_files[i]
  for(j in 1:nrow(mpileupout)){
    bases <- mpileupout[j,5]
    qnames <- mpileupout[j,6]
    dedup_res <- dedup_overlap(bases, qnames)
    mpileupout[j,5] <- paste(dedup_res[,2], collapse = "")
    mpileupout[j,6] <- paste(dedup_res[,1], collapse = ",")
    write.table(mpileupout, file = mpileupoutfile, sep = "\t", quote = F, row.names = F, col.names = F)
  }
}
