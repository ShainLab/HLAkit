#!/usr/bin/env Rscript
#!/opt/homebrew/bin Rscript
#!/opt/homebrew/share/man/man1 Rscript.1

# find HLA allele with genomic sequence using HLA-LA and HLA-HD result
library(optparse)

option_list <- list( 
    make_option(c("-l", "--hla_la"),
        type = "character",
        help = "HLA-LA result"
        ),
    make_option(c("-h", "--hla_hd"),
        type = "character",
        help = "HLA-HD result"
        ),
    make_option(c("-g", "--hladict"),
        type = "character",
        help = "Path to hla.dict file"
        ),
    make_option(c("-r", "--resultdir"),
        type = "character",
        help = "Output directory"
        )
)

# read the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# allot inputs to variables
hla_hd <- opt$hla_hd
hla_la <- opt$hla_la
genomic_alleles <- opt$hladict
resultdir <- opt$resultdir

if (is.null(hla_hd)) stop("ERROR: --hla_hd is required but not provided.")
if (is.null(hla_la)) stop("ERROR: --hla_la is required but not provided.")
if (is.null(genomic_alleles)) stop("ERROR: --hladict is required but not provided.")
if (is.null(resultdir)) stop("ERROR: --resultdir is required but not provided.")

set.seed(1000)
setwd(resultdir)

hd <- read.csv(hla_hd, sep = "\t", header = F)[,1]
la <- read.csv(hla_la, sep = "\t", header = F)[,1]
genomic <- read.csv(genomic_alleles, sep = '\t', header = F)

genomic <- genomic[-(1:2), 2]
genomic <- gsub('SN:', "", genomic)

# named vectors of hla alleles
la <- setNames(la, paste0("V", 1:6))
hd <- setNames(hd, paste0("V", 1:6))


match2digit <- function(allele) {
  parts <- strsplit(allele, "_")[[1]]
  paste(parts[1:3], collapse = "_")
}

hla_resolution <- function(allele) {
  length(strsplit(allele, "_")[[1]]) - 2
}

hla_truncate <- function(allele, res) {
  parts <- strsplit(allele, "_")[[1]]
  paste(parts[1:(2 + res)], collapse = "_")
}

find_genomic_allele <- function(la_allele, hd_allele, genomic) {
  la_base <- hla_resolution(la_allele)
  hd_base <- hla_resolution(hd_allele)
  
  for (res_offset in c(0, -1, -2)) {
    la_res <- la_base + res_offset
    hd_res <- hd_base + res_offset
    
    if (la_res > 0) {
      match <- grep(paste(hla_truncate(la_allele, la_res), "_", sep=""), genomic, fixed = TRUE)[1]
      if (!is.na(match)) return(genomic[match])
    }
    if (hd_res > 0) {
      match <- grep(paste(hla_truncate(hd_allele, hd_res), "_", sep=""), genomic, fixed = TRUE)[1]
      if (!is.na(match)) return(genomic[match])
    }
  }
  return(NA)
}

# match HLA-HD alleles to HLA-LA at upto 2-digit resolution
combos <- list(
  c(1,1,"a"), c(1,2,"a"), c(2,1,"a"), c(2,2,"a"),
  c(3,3,"b"), c(3,4,"b"), c(4,3,"b"), c(4,4,"b"),
  c(5,5,"c"), c(5,6,"c"), c(6,5,"c"), c(6,6,"c")
)

match_result <- lapply(combos, function(x) {
  la_idx <- as.integer(x[1])
  hd_idx <- as.integer(x[2])
  if (match2digit(la[la_idx]) == match2digit(hd[hd_idx])) hd[hd_idx] else "2_digit_match_not_found"
})

nr <- "2_digit_match_not_found"

# function to find allele with genomic sequence
resolve_gene <- function(m1, m2, m3, m4, raw1, raw2) {
  a1 <- NA; a2 <- NA

  if (m1 != nr) a1 <- m1
  if (m2 != nr) a1 <- m2
  if (m3 != nr) a2 <- m3
  if (m4 != nr) a2 <- m4
  if (m1 != nr && m2 != nr) { a1 <- m1; a2 <- m2 }
  if (m3 != nr && m4 != nr) { a1 <- m3; a2 <- m4 }

  if (is.na(a1) && is.na(a2)) {
    a1 <- raw1; a2 <- raw2
  } else if (is.na(a1) && !is.na(a2)) {
    if (m3 == nr) a1 <- raw1
    if (m4 == nr) a1 <- raw2
  } else if (!is.na(a1) && is.na(a2)) {
    if (m1 == nr) a2 <- raw1
    if (m2 == nr) a2 <- raw2
  }

  c(a1, a2)
}

hla_a_hd  <- resolve_gene(match_result[[1]], match_result[[2]], match_result[[3]], match_result[[4]], hd[1], hd[2])
hla_b_hd  <- resolve_gene(match_result[[5]], match_result[[6]], match_result[[7]], match_result[[8]], hd[3], hd[4])
hla_c_hd  <- resolve_gene(match_result[[9]], match_result[[10]], match_result[[11]], match_result[[12]], hd[5], hd[6])

la_new <- data.frame(
  gene  = c("A", "B", "C"),
  allele1 = c(la[1], la[3], la[5]),
  allele2 = c(la[2], la[4], la[6]),
  stringsAsFactors = FALSE
)

hd_new <- data.frame(
  gene    = c("A", "B", "C"),
  allele1 = c(hla_a_hd[1], hla_b_hd[1], hla_c_hd[1]),
  allele2 = c(hla_a_hd[2], hla_b_hd[2], hla_c_hd[2]),
  stringsAsFactors = FALSE
)

# add homozygous
fix_homozygous <- function(df) {
  for (j in 1:nrow(df)) {
    if (!is.na(df$allele1[j]) && df$allele1[j] == "homozygous") {
      df$allele1[j] <- df$allele2[j]
      df$allele2[j] <- "homozygous"
    }
    if (!is.na(df$allele1[j]) && !is.na(df$allele2[j]) && df$allele1[j] == df$allele2[j]) {
      df$allele2[j] <- "homozygous"
    }
  }
  df
}

la_new <- fix_homozygous(la_new)
hd_new <- fix_homozygous(hd_new)

genomic_alleles_new <- data.frame(
  gene    = c("A", "B", "C"),
  allele1 = NA_character_,
  allele2 = NA_character_,
  stringsAsFactors = FALSE
)

for (j in 1:3) {
  for (allele_col in 1:2) {
    la_al <- la_new[j, allele_col + 1]
    hd_al <- hd_new[j, allele_col + 1]

    if (!is.na(la_al) && la_al == "homozygous") {
      genomic_alleles_new[j, allele_col + 1] <- "homozygous"
    } else {
      hit <- find_genomic_allele(la_al, hd_al, genomic)
      genomic_alleles_new[j, allele_col + 1] <- if (!is.na(hit)) hit else "genomic_allele_not_found"
    }
  }
}

write.table(genomic_alleles_new, file = "genomic_alleles.txt", row.names = F, col.names = T, quote = F)
