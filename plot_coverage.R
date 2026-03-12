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
resultdir <- opt$resultdir
check_file <- function(f) {
    if (!file.exists(f)) stop(paste("File not found:", f))
    if (file.info(f)$size == 0) stop(paste("File is empty:", f))
}

check_file(coverage)
setwd(resultdir)

# read data
df <- read.delim(coverage, sep = "\t", header = T)

#normalize coverage
df$normal <- df$normal/as.numeric(WESnormalcoverage)
df$normal[which(is.nan(df$normal))] <- 0
df$tumor <- df$tumor/as.numeric(WEStumorcoverage)
df$tumor[which(is.nan(df$tumor))] <- 0

colnames(df)[2:3] <- c("Normal", "Tumor")

# Parse gene from allele
parse_allele <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  gene  <- tolower(parts[2])
  allele_num <- paste(parts[3:length(parts)], collapse = "_")
  list(gene = gene, label = allele_num, homozygous = FALSE)
}

parsed         <- lapply(df$allele, parse_allele)
df$gene        <- sapply(parsed, `[[`, "gene")
df$label       <- sapply(parsed, `[[`, "label")
df$homozygous  <- sapply(parsed, `[[`, "homozygous")

# Skip homozygous gene(s)
gene_counts      <- table(df$gene)
homozygous_genes <- names(gene_counts[gene_counts < 2])
df <- df[!df$gene %in% homozygous_genes, ]


if(nrow(df) == 0){
  stop("All genes are homozygous. Cannot plot allelic coverage.")
  } else{     
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

  # order by descending normal coverage
  normal_order <- df[order(df$gene, -df$Normal), c("gene", "label")]
  gene_allele_order <- split(normal_order$label, normal_order$gene)
  long$label <- as.character(long$label)

  pal <- c("Normal" = "#4E79A7", "Tumor" = "#E15759")

  # One plot per gene
  genes <- unique(long$gene)

  plot_list <- lapply(genes, function(g) {
    sub <- long[long$gene == g, ]
    alleles <- gene_allele_order[[g]]
    sub$label <- factor(sub$label, levels = alleles)
    allele_labels <- setNames(paste0("hla_", g, "_", alleles), alleles)
    names(allele_labels) <- alleles

    ggplot(sub, aes(x = label, y = coverage, fill = sample)) +
      geom_bar(
        stat     = "identity",
        position = position_dodge(width = 0.65),
        width    = 0.6,
        color    = "white",
        linewidth = 0.3
      ) +
      scale_fill_manual(values = pal, name = "Sample") +
      scale_x_discrete(
        labels = allele_labels
      ) +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.15)),
        limits = c(0, 1)
      ) +
      labs(
        title    = paste0("HLA-", toupper(g), " Coverage"),
        x        = NULL,
        y        = "Coverage"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title      = element_text(face = "bold", size = 14, hjust = 0),
        axis.text.x     = element_text(size = 8, angle = 20, hjust = 1),
        axis.text.y     = element_text(size = 10, color = "grey30"),
        axis.title.y    = element_text(size = 11, margin = margin(r = 8)),
        axis.line       = element_line(color = "grey70"),
        axis.ticks      = element_line(color = "grey70"),
        legend.position = "top",
        legend.title    = element_text(size = 10, face = "bold"),
        legend.text     = element_text(size = 10),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
        plot.margin     = margin(12, 16, 12, 12)
      )
  })

  pdf(paste0(resultdir, "/hla_coverage_plots.pdf"), width = 5, height = 4)
  invisible(lapply(plot_list, print))
  dev.off()

  message("Saved: hla_coverage_plots.pdf")

}

