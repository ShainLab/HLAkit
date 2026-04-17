#!/usr/bin/env Rscript
#!/opt/homebrew/bin Rscript
#!/opt/homebrew/share/man/man1 Rscript.1

# arm level loh
library(optparse)

option_list <- list( 
    make_option(c("-s", "--samplename"),
        type = "character",
        help = "Sample name"
        ),
    make_option(c("-l", "--loh"),
        type = "character",
        help = "Path to LOH file"
        ),
    make_option(c("-r", "--resultdir"),
        type = "character",
        help = "Result directory"),
    make_option(c("-p", "--tumor_purity"),
        type = "numeric",
        help = "Tumor purity in fraction")
)

# read the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# allot inputs to variables
samplename <- opt$samplename
loh <- opt$loh
resultdir <- opt$resultdir
tumor_purity <- opt$tumor_purity

if (is.null(samplename)) stop("ERROR: --samplename is required but not provided.")
if (is.null(loh)) stop("ERROR: --loh is required but not provided.")
if (is.null(resultdir)) stop("ERROR: --resultdir is required but not provided.")
if (is.null(tumor_purity)) stop("ERROR: --tumor_purity is required but not provided.")

suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DNAcopy))

setwd(resultdir)

merge_small_segments <- function(segments, df, vaf_col, chr_col, 
                                 pos_col, min_variants) {
  # Process each chromosome separately
  merged_list <- list()
  
  for(chrom in unique(segments$chr)) {
    chr_segs <- segments[segments$chr == chrom, ]
    chr_segs <- chr_segs[order(chr_segs$startpos), ]
    
    if(nrow(chr_segs) == 1) {
      merged_list[[length(merged_list) + 1]] <- chr_segs
      next
    }
    
    # Iteratively merge small segments
    repeat {
      small_idx <- which(chr_segs$numberofvariants < min_variants)
      
      if(length(small_idx) == 0) break
      
      # Process first small segment
      idx <- small_idx[1]
      
      # Determine which neighbor to merge with
      if(idx == 1) {
        # First segment - merge with next
        merge_with <- 2
      } else if(idx == nrow(chr_segs)) {
        # Last segment - merge with previous
        merge_with <- idx - 1
      } else {
        # Middle segment - merge with closer neighbor (by median VAF)
        prev_diff <- abs(chr_segs$median_vaf[idx] - chr_segs$median_vaf[idx-1])
        next_diff <- abs(chr_segs$median_vaf[idx] - chr_segs$median_vaf[idx+1])
        
        if(prev_diff < next_diff) {
          merge_with <- idx - 1
        } else {
          merge_with <- idx + 1
        }
      }
      
      # Merge segments
      new_startpos <- min(chr_segs$startpos[c(idx, merge_with)])
      new_endpos <- max(chr_segs$endpos[c(idx, merge_with)])
      
      # Get all variants in merged segment
      chr_name <- chr_segs$chr[idx]
      mask <- (df[[chr_col]] == paste0("chr", chr_name) | 
               df[[chr_col]] == chr_name) & 
              df[[pos_col]] >= new_startpos & 
              df[[pos_col]] <= new_endpos
      
      merged_vaf <- df[[vaf_col]][mask]
      
      # Create merged segment
      merged_seg <- data.frame(
        sample = chr_segs$sample[idx],
        chr = chr_segs$chr[idx],
        startpos = new_startpos,
        endpos = new_endpos,
        numberofvariants = length(merged_vaf),
        mean_vaf = mean(merged_vaf, na.rm = TRUE),
        median_vaf = median(merged_vaf, na.rm = TRUE)
      )
      
      # Remove old segments and add merged one
      keep_idx <- setdiff(1:nrow(chr_segs), c(idx, merge_with))
      chr_segs <- rbind(chr_segs[keep_idx, ], merged_seg)
      chr_segs <- chr_segs[order(chr_segs$startpos), ]
    }
    
    merged_list[[length(merged_list) + 1]] <- chr_segs
  }
  
  return(bind_rows(merged_list))
}


segment_with_cbs_merge <- function(df, vaf_col, chr_col, 
                                   pos_col, sample_id,
                                   min_variants, cluster_merge_threshold, thresh) {
  
  # Prepare data for DNAcopy
  cna_data <- data.frame(
    chrom = df[[chr_col]],
    maploc = df[[pos_col]],
    value = df[[vaf_col]]
  )
  
  # Remove chr prefix if present
  cna_data$chrom <- gsub("chr", "", cna_data$chrom)
  
  # Create CNA object
  CNA_object <- CNA(genomdat = cna_data$value,
                    chrom = cna_data$chrom,
                    maploc = cna_data$maploc,
                    data.type = "logratio",
                    sampleid = sample_id)
  
  # Smooth outliers
  CNA_smoothed <- smooth.CNA(CNA_object)
  
  # Segment
  segments_raw <- segment(CNA_smoothed, verbose = 1)
  
  # Extract segment information
  segments <- segments_raw$output
  colnames(segments) <- c("sample", "chr", "startpos", "endpos", 
                         "numberofvariants", "mean_vaf")
  
  segments <- segments[!is.na(segments$mean_vaf),,drop=F]
  # remove segments with less than or equal to 5 variants
  segments <- segments[segments$numberofvariants>5,]

  # Add median calculation
  segments$median_vaf <- NA
  for(i in 1:nrow(segments)) {
    chr_name <- segments$chr[i]
    # Handle both "chr6" and "6" formats
    mask <- (df[[chr_col]] == paste0("chr", chr_name) | df[[chr_col]] == chr_name) & 
            df[[pos_col]] >= segments$startpos[i] & 
            df[[pos_col]] <= segments$endpos[i]
    segments$median_vaf[i] <- median(df[[vaf_col]][mask], na.rm = TRUE)
  }
  
  
  # Merge small segments
  segments_merged <- merge_small_segments(segments, df, vaf_col, chr_col, 
                                         pos_col, min_variants)
  
  
  # Cluster segments by median VAF
  gmm_fit <- Mclust(segments_merged$median_vaf, G = 1:10)
  segments_merged$clusternum <- gmm_fit$classification
  
  cluster_medians <- tapply(segments_merged$median_vaf, 
                            segments_merged$clusternum, median)

  # Reorder cluster numbers
  cluster_medians <- sort(cluster_medians)
  # Create mapping: old cluster IDs -> new sequential IDs (1, 2, 3, ...)
  renamecluster <- setNames(1:length(cluster_medians), names(cluster_medians))

  # Apply mapping
  segments_merged$clusternum <- unname(renamecluster[as.character(segments_merged$clusternum)])


  # if cluster medians differ by less than threshold BAF, merge them.
   # Recalculate after reordering
  cluster_medians <- tapply(segments_merged$median_vaf, 
                            segments_merged$clusternum, median)
  cluster_medians <- sort(cluster_medians)
  
  
  # Merge clusters with close medians
  n_clusters <- length(cluster_medians)
  merge_map <- setNames(1:n_clusters, names(cluster_medians))
  
  if(n_clusters > 1) {
    new_cluster_id <- 1
    merge_map[1] <- new_cluster_id
    
    for(i in 2:n_clusters) {
      diff_from_prev <- cluster_medians[i] - cluster_medians[i-1]
      
      if(diff_from_prev < cluster_merge_threshold) {
        merge_map[i] <- merge_map[i-1]
        
      } else {
        new_cluster_id <- new_cluster_id + 1
        merge_map[i] <- new_cluster_id
      }
    }
    
    # Apply merge
    old_clusters <- segments_merged$clusternum
    segments_merged$clusternum <- merge_map[as.character(old_clusters)]
  }
  
  # Final cluster medians
  final_cluster_medians <- tapply(segments_merged$median_vaf, 
                                  segments_merged$clusternum, median)
  final_cluster_medians <- sort(final_cluster_medians)

  segments <- segments_merged
  segments <- segments[!is.na(segments$clusternum),,drop=F]
  if(all(segments$clusternum > 1)) segments$clusternum <- segments$clusternum-1

  # assigning clusternum based on closeness to medians instead of mclust
  for(n in 1:nrow(segments)){
    segments$clusternum[n] <- which.min(abs(final_cluster_medians - segments$median_vaf[n]))
  }

  hla_a <- max(segments$clusternum[which(segments$chr==6 & segments$startpos <= 29913661 & segments$endpos >= 29909037)])
  hla_b <- max(segments$clusternum[which(segments$chr==6 & segments$startpos <= 31324964 & segments$endpos >= 31321649)])
  hla_c <- max(segments$clusternum[which(segments$chr==6 & segments$startpos <= 31239869 & segments$endpos >= 31236526)])

  if(hla_a>1) out <- "TRUE" else out <- "FALSE"
  if(hla_b>1) out <- c(out,"TRUE") else out <- c(out,  "FALSE")
  if(hla_c>1) out <- c(out,"TRUE") else out <- c(out,  "FALSE")

  names(out) <- c("hla_a", "hla_b", "hla_c")

  signal <-  which.min(abs(final_cluster_medians - thresh))
  fcm <- unname(final_cluster_medians)
  
  if(fcm[signal] > fcm[1]){
    noisy <- FALSE
  } else if(fcm[signal] > fcm[1] && length(fcm) == 1){
    noisy <- FALSE
  } else{
    noisy <- TRUE
  }
  names(noisy) <- "Noisy_data"
  return(list(out = out, noisy = noisy, segments = segments, gmm = gmm_fit, cbs = segments_raw, cluster_medians = final_cluster_medians))
}


outfile <- paste(samplename, ".hla_loh_segments.txt", sep = "")
noisyoutfile <- sub(".txt", ".noisy.txt", outfile)

# VAF data
df <- read.delim(loh, sep = "\t")
df <- unique(df)
df$Chromosome <- gsub("chr", "", df$Chromosome)
vaf_col <- "Tumor_VAF_corrected"
chr <- "Chromosome"
pos <- "Start_Position"
index <- which(df$Chromosome=="X")
if(length(index)>0) df <- df[-index,]

cat("Sample:", samplename, "\n")
thresh <- tumor_purity/2
cat("Finding segments and clusters \n")
result <- segment_with_cbs_merge(df, vaf_col = vaf_col, chr_col = chr, 
                                 pos_col = pos, sample_id=samplename, min_variants = 30, cluster_merge_threshold = 0.053, thresh = thresh)

write.table(result$out, outfile, sep = "\t", quote=F, col.names = F)
write.table(result$noisy, noisyoutfile, sep = "\t", quote=F, col.names = F, row.names = F)
write.table(result$segments, sub(".txt", ".seg", outfile), row.names = F, sep = "\t", quote=F, col.names = F)
