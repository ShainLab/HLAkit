#!/usr/bin/env Rscript
#!/opt/homebrew/bin Rscript
#!/opt/homebrew/share/man/man1 Rscript.1

# hla region arm level copy number alteration

library(mclust)
library(dplyr)
library(optparse)
library(ggplot2)

option_list <- list( 
    make_option(c("-s", "--samplename"),
        type = "character",
        help = "Sample name"
        ),
    make_option(c("-c", "--cns"),
        type = "character",
        help = "Path to CNS file"
        ),
    make_option(c("-c", "--cnr"),
        type = "character",
        help = "Path to CNR file"
        ),
    make_option(c("-r", "--resultdir"),
        type = "character",
        help = "Result directory"),
    make_option(c("-r", "--tumor_purity"),
        type = "numeric",
        help = "Tumor purity in fraction")
)

# read the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# allot inputs to variables
samplename <- opt$samplename
cns <- opt$cns
cnr <- opt$cnr
resultdir <- opt$resultdir
tumor_purity <- opt$tumor_purity

if (is.null(samplename)) stop("ERROR: --samplename is required but not provided.")
if (is.null(cns)) stop("ERROR: --cns is required but not provided.")
if (is.null(cnr)) stop("ERROR: --cnr is required but not provided.")
if (is.null(resultdir)) stop("ERROR: --resultdir is required but not provided.")
if (is.null(tumor_purity)) stop("ERROR: --tumor_purity is required but not provided.")



expected_ratios <- function(cn, purity, ploidy) {
  return(log2((purity * cn + (1 - purity) * ploidy) / ploidy))
}

cnv_clusters <- function(segments, tumor_purity, ploidy = 2) {
  
  # Calculate expected log2 ratios for different CN states
  # Formula: log2((purity * CN + (1-purity) * ploidy) / ploidy)
  
  # Expected values for common copy number states
  expected_cn1 <- expected_ratios(1, tumor_purity, ploidy)  # Heterozygous loss
  expected_cn2 <- expected_ratios(2, tumor_purity, ploidy)  # Neutral
  expected_cn3 <- expected_ratios(3, tumor_purity, ploidy)  # Single copy gain
  
  # Fit GMM
  gmm <- Mclust(segments$log2, G = 3)
  
  cluster_means <- gmm$parameters$mean
  cluster_sds <- sqrt(gmm$parameters$variance$sigmasq)
  
  # Identify clusters based on proximity to expected values
  neutral_idx <- which.min(abs(cluster_means - expected_cn2))
  loss_idx <- which.min(cluster_means)
  gain_idx <- which.max(cluster_means)
  
  neutral_mean <- cluster_means[neutral_idx]
  neutral_sd <- cluster_sds[neutral_idx]
  
  # Set thresholds based on expected cn
  # Loss threshold: midpoint between expected neutral and expected loss
  loss_threshold_bio <- (expected_cn2 + expected_cn1) / 2
  loss_threshold <- loss_threshold_bio
  
  # use whichever is more conservative
  loss_threshold_data <- (neutral_mean + cluster_means[loss_idx]) / 2
  loss_threshold <- min(loss_threshold_bio, loss_threshold_data)
  
  # Gain threshold
  gain_threshold_bio <- (expected_cn2 + expected_cn3) / 2
  gain_threshold <- gain_threshold_bio

  # In very low purity samples, might need to relax slightly
  if (tumor_purity < 0.3) {
    cat("WARNING: Low tumor purity (<30%). Thresholds may be less reliable.\n")
    # For low purity, conservative threshold for gain
    gain_threshold <- max(gain_threshold, 0.15)
  }

  # Classify segments
  segments$cn_state <- case_when(
    segments$log2 < loss_threshold ~ "Loss_hi",
    segments$log2 < loss_threshold - 0.40*loss_threshold ~ "Loss_low",
    segments$log2 > gain_threshold ~ "Gain_hi",
    segments$log2 > gain_threshold - 0.40*gain_threshold ~ "Gain_low",
    TRUE ~ "Neutral"
  )
  # Add expected CN based on classification
  segments$estimated_cn <- case_when(
    segments$cn_state == "Loss_low" ~ 0.5,
    segments$cn_state == "Loss_hi" ~ 1,
    segments$cn_state == "Neutral" ~ 2,
    segments$cn_state == "Gain_low" ~ 3,
    segments$cn_state == "Gain_hi" ~ 3.5

  )
  
  hla_a <- segments$estimated_cn[which(segments$chr==6 & segments$startpos <= 29913661 & segments$endpos >= 29909037)]
  hla_b <- segments$estimated_cn[which(segments$chr==6 & segments$startpos <= 31324964 & segments$endpos >= 31321649)]
  hla_c <- segments$estimated_cn[which(segments$chr==6 & segments$startpos <= 31239869 & segments$endpos >= 31236526)]

  if(length(hla_a)>1){
    if(any(hla_a>=2)) hla_a <- "Gain" else if(any(hla_a<=2)) hla_a <- "Loss" else hla_a <- "Neutral"
    } else{
      if(hla_a==3.5) hla_a <- "Gain_hi" else if (hla_a==3) hla_a <- "Gain_low" else if(hla_a==1) hla_a <- "Loss_hi" else if(hla_a==0.5) hla_a <- "Loss_low" else hla_a <- "Neutral"
    }
  if(length(hla_b)>1){
    if(any(hla_b>=2)) hla_b <- "Gain" else if(any(hla_b<=2)) hla_b <- "Loss" else hla_b <- "Neutral"
    } else{
      if(hla_b==3.5) hla_b <- "Gain_hi" else if (hla_b==3) hla_b <- "Gain_low" else if(hla_b==1) hla_b <- "Loss_hi" else if(hla_b==0.5) hla_b <- "Loss_low" else hla_b <- "Neutral"
    }
  if(length(hla_c)>1){
    if(any(hla_c>=2)) hla_c <- "Gain" else if(any(hla_c<=2)) hla_c <- "Loss" else hla_c <- "Neutral"
    } else{
      if(hla_c==3.5) hla_c <- "Gain_hi" else if (hla_c==3) hla_c <- "Gain_low" else if(hla_c==1) hla_c <- "Loss_hi" else if(hla_c==0.5) hla_c <- "Loss_low" else hla_c <- "Neutral"
    }
  out <- c(hla_a, hla_b, hla_c)
  names(out) <- c("hla_a", "hla_b", "hla_c")

  return(list(
    segments = segments,
    out = out,
    thresholds = c(loss = loss_threshold, loss_low = loss_threshold - 0.40*loss_threshold, gain = gain_threshold, gain_low = gain_threshold - 0.40*gain_threshold),
    expected_values = c(cn1 = expected_cn1, 
                       cn2 = expected_cn2, cn3 = expected_cn3),
    model = gmm
  ))
}


setwd(resultdir)

# run cnv cluster
cns <- read.delim(cns)
cnr <- read.delim(cnr)
colnames(cns) <- c("chr", "startpos", "endpos", "gene", "log2", "depth", "numberofvariants", "weight")
cns <- cns[,-4]
cns$chr <- gsub("chr", "", cns$chr)
cns <- cns[which(cns$chr != "X" & cns$chr != "Y"),]
cns <- as.data.frame(cbind("sample"=samplename, cns))

#delete segments with numberofvariants less than 30
cns <- cns[which(cns$numberofvariants > 30),]
result <- cnv_clusters(cns, tumor_purity=tumor_purity)

outname <- paste(samplename, ".cns_clusters.txt", sep ="")
write.table(result$segments, file=outname, sep ='\t', quote=F, row.names=F )
hlaoutname <- paste(samplename, ".hla.cn.txt", sep ="")
write.table(result$out, file=hlaoutname, sep ='\t', quote=F, col.names = F )

