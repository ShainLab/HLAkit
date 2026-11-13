#!/usr/bin/env Rscript
# find HLA class I allele with genomic sequence using HLA-LA and HLA-HD result

library(optparse)
library(openxlsx)

option_list <- list( 
    make_option(c("-l", "--hla_la"), type = "character", help = "HLA-LA result"),
    make_option(c("-d", "--hla_hd"), type = "character", help = "HLA-HD result"),
    make_option(c("-g", "--hladict"), type = "character", help = "Path to hla.dict file"),
    make_option(c("-r", "--resultdir"), type = "character", help = "Output directory"),
    make_option(c("-o", "--outfile"), type = "character", help = "Output xlsx filename")
)

opt <- parse_args(OptionParser(option_list = option_list))

hla_hd <- opt$hla_hd
hla_la <- opt$hla_la
genomic_alleles <- opt$hladict
resultdir <- opt$resultdir
outfile <- opt$outfile

if (is.null(hla_hd)) stop("ERROR: --hla_hd is required but not provided.")
if (is.null(hla_la)) stop("ERROR: --hla_la is required but not provided.")
if (is.null(genomic_alleles)) stop("ERROR: --hladict is required but not provided.")
if (is.null(resultdir)) stop("ERROR: --resultdir is required but not provided.")
if (is.null(outfile)) stop("ERROR: --outfile is required but not provided.")

# hla_la file should contain hla_la results for all samples without col names; rownames=sample
# hla_hd file should contain hla_hd results for all samples without col names; rownames=sample

set.seed(1000)

setwd(resultdir)

hla_la_all <- read.csv(hla_la, sep = "\t", header = F)
hla_hd_all <- read.csv(hla_hd, sep = "\t", header = F)
genomic <- read.csv(genomic_alleles, sep = '\t', header = FALSE)
genomic <- genomic[-(1:2), 2]
genomic <- gsub('SN:', "", genomic)

genes <- c("hla_a", "hla_b", "hla_c")

# ---- gene detection (order-independent) ----
extract_gene <- function(allele) {
  if (is.na(allele) || allele == "homozygous") return(NA_character_)
  parts <- strsplit(allele, "_")[[1]]
  paste(parts[1:2], collapse = "_")
}

# group a raw allele vector into a named list of up to-2-allele vectors per gene
group_by_gene <- function(raw, genes) {
  raw_genes <- sapply(raw, extract_gene)
  setNames(lapply(genes, function(g) {
    hits <- raw[which(raw_genes == g)]
    length(hits) <- 2   # pad with NA if a gene has <2 calls
    hits
  }), genes)
}


match2digit <- function(allele) {
  if (is.na(allele)) return(NA_character_)
  parts <- strsplit(allele, "_")[[1]]
  if (length(parts) < 3) return(NA_character_)
  paste(parts[1:3], collapse = "_")
}

hla_resolution <- function(allele) {
  if (is.na(allele)) return(0)
  length(strsplit(allele, "_")[[1]]) - 2
}

hla_truncate <- function(allele, res) {
  parts <- strsplit(allele, "_")[[1]]
  paste(parts[1:(2 + res)], collapse = "_")
}

find_genomic_allele <- function(la_allele, hd_allele, genomic) {
  if (is.na(la_allele) && is.na(hd_allele)) return(NA)
  la_base <- hla_resolution(la_allele)
  hd_base <- hla_resolution(hd_allele)

  for (res_offset in c(0, -1, -2)) {
    la_res <- la_base + res_offset
    hd_res <- hd_base + res_offset

    if (!is.na(la_allele) && la_res > 0) {
      match <- grep(paste(hla_truncate(la_allele, la_res), "_", sep = ""), genomic, fixed = TRUE)[1]
      if (!is.na(match)) return(genomic[match])
    }
    if (!is.na(hd_allele) && hd_res > 0) {
      match <- grep(paste(hla_truncate(hd_allele, hd_res), "_", sep = ""), genomic, fixed = TRUE)[1]
      if (!is.na(match)) return(genomic[match])
    }
  }
  return(NA)
}


nr <- "2_digit_match_not_found"

# 2-digit match across all 4 la/hd allele combos, per gene
match_gene <- function(la_pair, hd_pair) {
  combos <- list(c(1, 1), c(1, 2), c(2, 1), c(2, 2))
  lapply(combos, function(x) {
    la_al <- la_pair[x[1]]
    hd_al <- hd_pair[x[2]]
    if (!is.na(la_al) && !is.na(hd_al) && match2digit(la_al) == match2digit(hd_al)) {
      hd_al
    } else {
      nr
    }
  })
}

resolve_gene <- function(m1, m2, m3, m4, raw1, raw2) {
  a1 <- NA; a2 <- NA
  if (m1 != nr) a1 <- m1
  if (m2 != nr) a1 <- m2
  if (m3 != nr) a2 <- m3
  if (m4 != nr) a2 <- m4
  if (m1 != nr && m2 != nr) { a1 <- m1; a2 <- m2 }
  if (m3 != nr && m4 != nr) { a1 <- m3; a2 <- m4 }

  if (!is.na(a1) && !is.na(a2) && a1 == a2 && raw1 != raw2) {
    a1 <- raw1; a2 <- raw2
  }

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

fix_homozygous_raw <- function(raw) {
  for (i in seq(1, length(raw), by = 2)) {
    a1 <- raw[i]; a2 <- raw[i + 1]
    if (!is.na(a1) && a1 == "homozygous" && !is.na(a2)) raw[i] <- a2
    if (!is.na(a2) && a2 == "homozygous" && !is.na(a1)) raw[i + 1] <- a1
  }
  raw
}

result <- lapply(1:nrow(hla_la_all), function(z) {
  la_raw <- unlist(hla_la_all[z, 2:7], use.names = FALSE)
  hd_raw <- unlist(hla_hd_all[z, 2:7], use.names = FALSE)

  la_raw <- fix_homozygous_raw(la_raw)
  hd_raw <- fix_homozygous_raw(hd_raw)

  # group this sample's alleles by gene (order-independent)
  la_by_gene <- group_by_gene(la_raw, genes)
  hd_by_gene <- group_by_gene(hd_raw, genes)

  hd_resolved <- lapply(genes, function(g) {
    la_pair <- la_by_gene[[g]]
    hd_pair <- hd_by_gene[[g]]
    m <- match_gene(la_pair, hd_pair)
    resolve_gene(m[[1]], m[[2]], m[[3]], m[[4]], hd_pair[1], hd_pair[2])
  })
  names(hd_resolved) <- genes

  la_new <- data.frame(
    gene    = genes,
    allele1 = sapply(genes, function(g) la_by_gene[[g]][1]),
    allele2 = sapply(genes, function(g) la_by_gene[[g]][2]),
    stringsAsFactors = FALSE
  )

  hd_new <- data.frame(
    gene    = genes,
    allele1 = sapply(genes, function(g) hd_resolved[[g]][1]),
    allele2 = sapply(genes, function(g) hd_resolved[[g]][2]),
    stringsAsFactors = FALSE
  )

  la_new <- fix_homozygous(la_new)
  hd_new <- fix_homozygous(hd_new)

  genomic_alleles_new <- data.frame(
    gene    = genes,
    allele1 = NA_character_,
    allele2 = NA_character_,
    stringsAsFactors = FALSE
  )

  for (j in seq_along(genes)) {
    for (allele_col in 1:2) {
      la_al <- la_new[j, allele_col + 1]
      hd_al <- hd_new[j, allele_col + 1]

      if ((!is.na(la_al) && la_al == "homozygous") || (!is.na(hd_al) && hd_al == "homozygous")) {
        genomic_alleles_new[j, allele_col + 1] <- "homozygous"
      } else {
        hit <- find_genomic_allele(la_al, hd_al, genomic)
        genomic_alleles_new[j, allele_col + 1] <- if (!is.na(hit)) hit else "genomic_allele_not_found"
      }
    }
  }
  return(genomic_alleles_new)
})

# flatten each sample's 3-gene x 2-allele result into one row
result1 <- lapply(1:length(result), function(x) {
  unlist(lapply(1:length(genes), function(j) unlist(result[[x]][j, 2:3], use.names = FALSE)))
})

result2 <- do.call('rbind', result1)

col_names <- unlist(lapply(genes, function(g) c(paste0(g, "_allele1"), paste0(g, "_allele2"))))
colnames(result2) <- col_names
rownames(result2) <- hla_la_all[,1]

write.xlsx(result2, file = outfile, rowNames=T)
