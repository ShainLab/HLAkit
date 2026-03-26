#!/usr/bin/env Rscript
#!/opt/homebrew/bin Rscript
#!/opt/homebrew/share/man/man1 Rscript.1

# Find loss of heterozygosity

library(optparse)
library(ggplot2)

option_list <- list( 
    make_option(c("-c", "--coverage"),
        type = "character",
        help = "Normal and tumor coverage of HLA alleles"
        ),
    make_option(c("-n", "--WESnormalcoverage"),
        type = "numeric",
        help = "WES normal coverage"
        ),
    make_option(c("-t", "--WEStumorcoverage"),
        type = "numeric",
        help = "WES tumor coverage"
        ),
    make_option(c("--aibfile"),
        type = "character",
        help = "Allelic imbalance limits file"
        ),
    make_option(c("--lowcov"),
        type = "character",
        help = "Low coverage limits file"
        ),
    make_option(c("-r", "--resultdir"),
        type = "character",
        help = "Result directory")
)
# read the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# allot inputs to variables
coverage <- opt$coverage
WESnormalcoverage <- opt$WESnormalcoverage
WEStumorcoverage <- opt$WEStumorcoverage
aibfile <- opt$aibfile
lowcov <- opt$lowcov
resultdir <- opt$resultdir


if (is.null(coverage)) stop("ERROR: --coverage is required but not provided.")
if (is.null(WESnormalcoverage)) stop("ERROR: --WESnormalcoverage is required but not provided.")
if (is.null(WEStumorcoverage)) stop("ERROR: --WEStumorcoverage is required but not provided.")
if (is.null(aibfile)) stop("ERROR: --aibfiles is required but not provided.")
if (is.null(lowcov)) stop("ERROR: --lowcov is required but not provided.")
if (is.null(resultdir)) stop("ERROR: --resultdir is required but not provided.")

if (WESnormalcoverage == 0) stop("ERROR: --WESnormalcoverage is 0; division by zero would produce NaN for all normal coverage values.")
if (WEStumorcoverage  == 0) stop("ERROR: --WEStumorcoverage is 0; division by zero would produce NaN for all tumor coverage values.")

if (!dir.exists(resultdir)) stop(paste("ERROR: Result directory does not exist:", resultdir))

check_file <- function(f) {
    if (!file.exists(f)) stop(paste("File not found:", f))
    if (file.info(f)$size == 0) stop(paste("File is empty:", f))
}

check_file(coverage)
check_file(aibfile)
check_file(lowcov)
setwd(resultdir)

# read data
df <- read.delim(coverage, sep = "\t", header = T)
aib <- read.delim(aibfile, sep = "\t", header = T)
lowcov <- read.delim(lowcov, sep = "\t", header = T)

# file structure checks
if (nrow(df) == 0) stop(paste("ERROR: Coverage file has no data rows:", coverage))
if (nrow(aib) == 0) stop(paste("ERROR: AIB limits file has no data rows:", opt$aib))
if (nrow(lowcov) == 0) stop(paste("ERROR: Low coverage limits file has no data rows:", opt$lowcov))

required_cols_df <- c("allele", "normal", "tumor")
missing_df <- setdiff(required_cols_df, colnames(df))
if (length(missing_df) > 0) stop(paste("ERROR: Coverage file is missing required column(s):", paste(missing_df, collapse = ", ")))

if (nrow(aib) != nrow(lowcov)) {
  stop(paste0("ERROR: Row count mismatch — aib file has ", nrow(aib),
              " rows but lowcov limits file has ", nrow(lowcov), " rows. These must match."))
}

#normalize coverage
df$normal <- df$normal/as.numeric(WESnormalcoverage)
df$normal[which(is.nan(df$normal))] <- 0
df$tumor <- df$tumor/as.numeric(WEStumorcoverage)
df$tumor[which(is.nan(df$tumor))] <- 0

colnames(df)[2:3] <- c("Normal", "Tumor")

# Parse gene from allele
parse_allele <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  if (length(parts) < 3) {
    stop(paste0("ERROR: Allele name '", x, "' does not match expected format. ",
                "Got ", length(parts), " underscore-delimited part(s); expected at least 3."))
  }
  gene  <- tolower(parts[2])
  allele_num <- paste(parts[3:length(parts)], collapse = "_")
  list(gene = gene, label = allele_num, homozygous = FALSE)
}

parsed <- lapply(df$allele, parse_allele)
df$gene <- sapply(parsed, `[[`, "gene")
df$label <- sapply(parsed, `[[`, "label")
df$homozygous <- sapply(parsed, `[[`, "homozygous")

# Skip homozygous gene(s)
gene_counts      <- table(df$gene)
homozygous_genes <- names(gene_counts[gene_counts < 2])

if (length(homozygous_genes) > 0) {
  message("WARNING: The following gene(s) appear homozygous (only 1 allele detected) and will be excluded from plotting: ",
          paste(toupper(homozygous_genes), collapse = ", "))
}

homozygous_index <- which(df$gene %in% homozygous_genes)
if (length(homozygous_index) > 0) {
  aib    <- aib[-homozygous_index, 2:4]
  lowcov <- lowcov[-homozygous_index, 2:4]
} else {
  aib    <- aib[, 2:4]
  lowcov <- lowcov[, 2:4]
}

df <- df[!df$gene %in% homozygous_genes, ]

if(nrow(df) == 0){
  stop("ERROR: All genes are homozygous. Cannot plot allelic coverage.")
}

# aib/lowcov column structure after subsetting
if (ncol(aib) < 3)    stop(paste0("ERROR: AIB limits file has fewer than 3 columns after subsetting (found ", ncol(aib), "). Expected: flag, threshold, status."))
if (ncol(lowcov) < 3) stop(paste0("ERROR: Low coverage limits file has fewer than 3 columns after subsetting (found ", ncol(lowcov), "). Expected: flag, threshold, status."))

df$aibthresh <- as.numeric(aib[,2])
df$lowcovthresh <- as.numeric(lowcov[,2])

if (all(is.na(df$aibthresh))) {
  warning("WARNING: All AIB threshold values are NA after coercion. Check that column 2 of the AIB limits file contains numeric values.")
}
if (all(is.na(df$lowcovthresh))) {
  warning("WARNING: All low coverage threshold values are NA after coercion. Check that column 2 of the low coverage limits file contains numeric values.")
}

df$aibthreshval <- ifelse(aib[,1] == "Flag", "Flag",
  ifelse(aib[,3] == "Normal", "No",
  ifelse(aib[,3] == "Low", "Yes", "")))
df$lowcovthreshval <- ifelse(lowcov[,1] == "Flag", "Flag",
  ifelse(lowcov[,3] == "Normal", "No",
  ifelse(lowcov[,3] == "Low", "Yes", "")))
 
# reshape
long <- reshape(
  df[, c("label", "gene", "Normal", "Tumor")],
  varying   = c("Normal", "Tumor"),
  v.names   = "coverage",
  timevar   = "sample",
  times     = c("Normal", "Tumor"),
  direction = "long"
)
rownames(long) <- NULL

long$sample <- factor(long$sample, levels = c("Normal", "Tumor"))

tumor_sub <- long[long$sample == "Tumor", ]
tumor_sub <- as.data.frame(cbind(tumor_sub, df[, c("label", "aibthresh", "lowcovthresh", 
                                      "aibthreshval", "lowcovthreshval")]))

# Build threshold annotation data
thresh_long <- data.frame(
  label     = rep(tumor_sub$label, 2),
  gene      = rep(tumor_sub$gene, 2),
  threshold = rep(c("AIB", "LowCov"), each = nrow(tumor_sub)),
  value     = c(tumor_sub$aibthresh, tumor_sub$lowcovthresh),
  met       = c(tumor_sub$aibthreshval, tumor_sub$lowcovthreshval),
  row       = rep(c(1, 2), each = nrow(tumor_sub)),
  stringsAsFactors = FALSE
)

# order by descending normal coverage
gene_allele_order <- split(df$label, df$gene)
long$label <- as.character(long$label)

pal <- c("Normal" = "#4E79A7", "Tumor" = "#E15759")

# One plot per gene
genes <- unique(long$gene)

plot_list <- lapply(genes, function(g) {
  sub     <- long[long$gene == g, ]
  alleles <- gene_allele_order[[g]]
  sub$label <- factor(sub$label, levels = alleles)
 
  allele_labels <- setNames(paste0("hla_", g, "_", alleles), alleles)
 
  tg <- thresh_long[thresh_long$gene == g, ]
  tg$label <- factor(tg$label, levels = alleles)
 
  # y positions below x-axis, scaled per gene
  tg$y <- 0.15+max(sub$coverage) + (tg$row * 0.04)
  tg$yy <- max(tg$y) + 0.05

  #remove aib threshold from major allele
  tumor_cov <- sub$coverage[sub$sample == "Tumor"]
  if (length(tumor_cov) == 0 || all(is.na(tumor_cov))) {
    warning(paste0("WARNING: No tumor coverage values found for gene '", toupper(g), "'. Skipping AIB major-allele removal step for this gene."))
  } else {
    major_allele <- sub$label[sub$sample == "Tumor"][which.max(sub$coverage[sub$sample == "Tumor"])]
    index <- which(tg$threshold == "AIB" & tg$label == major_allele)
    if (length(index) == 0) {
      warning(paste0("WARNING: Could not find AIB threshold row for major allele '", major_allele,
                     "' in gene '", toupper(g), "'. The AIB annotation for this allele may be incorrect."))
    } else {
      tg <- tg[-index,]
    }
  }

  # change lowcov thresh to Flag for Flagged alleles
  tg$value <- round(tg$value, 2)
  if(all(tg$met == "Flag")) tg$value <- "Flag"

  ggplot(sub, aes(x = label, y = coverage, fill = sample)) +
    geom_bar(
      stat      = "identity",
      position  = position_dodge(width = 0.65),
      width     = 0.6,
      color     = "white",
      linewidth = 0.3
    ) +
    geom_text(
      aes(label = sprintf("%.2f", coverage)),
      position    = position_dodge(width = 0.65),
      vjust       = -0.5,
      size        = 3,
      color       = "grey30",
      inherit.aes = TRUE
    ) +
    geom_text(data = tg,
      aes(x = label, y = yy, label = "Threshold:"),
      size        = 2.5,
      hjust       = 0.5,
      inherit.aes = FALSE
      )+
    geom_text(
      data = tg,
      aes(x     = label,
          y     = y,
          label = paste0(threshold, "=", value, " (Met: ", met, ")")),
      size        = 2.5,
      hjust       = 0.5,
      inherit.aes = FALSE
    ) +
    scale_color_identity() +
    scale_fill_manual(values = pal, name = NULL) +
    scale_x_discrete(labels = allele_labels) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.15))
    ) +
    coord_cartesian(ylim = c(0, 0.25+max(sub$coverage)), clip = "off") +
    labs(
      title    = paste0("HLA-", toupper(g), " Coverage"),
      x        = NULL,
      y        = "Normalized Coverage (HLA/Exome-wide)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 14, hjust = 0),
      axis.text.x     = element_text(size = 8, angle = 20, hjust = 1),
      axis.text.y     = element_text(size = 10, color = "grey30"),
      axis.title.y    = element_text(size = 11, margin = margin(r = 8)),
      axis.line       = element_line(color = "grey70"),
      axis.ticks      = element_line(color = "grey70"),
      legend.position = "right",
      legend.title    = element_text(size = 10, face = "bold"),
      legend.text     = element_text(size = 10),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
      plot.margin = margin(12, 16, 40, 12)
    )
})

pdf(paste0(resultdir, "/hla_coverage_plots.pdf"), width = 5, height = 4)
invisible(lapply(plot_list, print))
dev.off()

message("Saved: hla_coverage_plots.pdf")

