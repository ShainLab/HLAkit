#!/usr/bin/env Rscript
#!/opt/homebrew/bin Rscript
#!/opt/homebrew/share/man/man1 Rscript.1

# Find loss of heterozygosity

library(optparse)
library(dplyr)

option_list <- list( 
    make_option(c("-p", "--tumor_purity"), 
        type = "numeric",
        help = "tumor purity in fraction"
        ),
    make_option(c("-n", "--WESnormalcoverage"), 
        type = "numeric",
        help = "Exome-wide coverage of normal sample"
        ),
    make_option(c("-t", "--WEStumorcoverage"), 
        type = "numeric",
        help = "Exome-wide coverage of tumor sample"
    		),
    make_option(c("-c", "--coverage"),
        type = "character",
        help = "Normal and tumor coverage of HLA alleles"
        ),
    make_option(c("-r", "--resultdir"),
        type = "character",
        help = "Result directory")
)
# read the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# allot inputs to variables
tumor_purity <- opt$tumor_purity
WESnormalcoverage <- opt$WESnormalcoverage
WEStumorcoverage <- opt$WEStumorcoverage
coverage <- opt$coverage
resultdir <- opt$resultdir

# --- Argument validation ---
if (is.null(tumor_purity)) stop("ERROR: -p/--tumor_purity is required but not provided.")
if (is.null(WESnormalcoverage)) stop("ERROR: -n/--WESnormalcoverage is required but not provided.")
if (is.null(WEStumorcoverage)) stop("ERROR: -t/--WEStumorcoverage is required but not provided.")
if (is.null(coverage)) stop("ERROR: -c/--coverage is required but not provided.")
if (is.null(resultdir)) stop("ERROR: -r/--resultdir is required but not provided.")

if (tumor_purity < 0 || tumor_purity > 1) stop(paste0("ERROR: --tumor_purity must be between 0 and 1 (got ", tumor_purity, "). Provide purity as a fraction, not a percentage."))
if (WESnormalcoverage <= 0) stop(paste0("ERROR: --WESnormalcoverage must be > 0 (got ", WESnormalcoverage, ")."))
if (WEStumorcoverage  <= 0) stop(paste0("ERROR: --WEStumorcoverage must be > 0 (got ", WEStumorcoverage, ")."))
if (!dir.exists(resultdir)) stop(paste("ERROR: Result directory does not exist:", resultdir))

check_file <- function(f) {
    if (!file.exists(f)) stop(paste("File not found:", f))
    if (file.info(f)$size == 0) stop(paste("File is empty:", f))
}

check_file(coverage)
setwd(resultdir)
outname <- paste0(resultdir, "/loh.txt")

data <- read.delim(coverage, sep = "\t", header = T)

# --- Post-read structural checks ---
if (nrow(data) == 0) stop(paste("ERROR: Coverage file has no data rows:", coverage))

required_cols <- c("allele", "normal", "tumor")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) stop(paste("ERROR: Coverage file is missing required column(s):", paste(missing_cols, collapse = ", ")))

# create vectors of allele, normal_cov and tumor_cov
data <- data %>%
  mutate(gene = sub("(hla_[abc])_.*", "\\1", allele))

# --- Check gene parsing worked ---
if (all(data$gene == data$allele)) {
  warning("WARNING: Gene parsing via '(hla_[abc])_.*' matched nothing — all 'gene' values equal the full allele name. Check that allele names follow the expected hla_a/b/c_... format.")
}

gene_counts <- table(data$gene)

# --- Warn about unexpected gene count ---
expected_genes <- c("hla_a", "hla_b", "hla_c")
found_genes <- names(gene_counts)
unexpected_genes <- setdiff(found_genes, expected_genes)
if (length(unexpected_genes) > 0) {
  warning(paste("WARNING: Unexpected gene names found after parsing:", paste(unexpected_genes, collapse = ", "),
                "— expected only hla_a, hla_b, hla_c. Check allele name format."))
}
homozygous_genes <- names(gene_counts[gene_counts < 2])
if (length(homozygous_genes) > 0) {
  message("NOTE: Homozygous gene(s) detected (only 1 allele in coverage file); placeholder row will be added: ",
          paste(toupper(homozygous_genes), collapse = ", "))
}

alleles <- character()
normal_cov <- numeric()
tumor_cov <- numeric()

for (g in names(gene_counts)) {
  rows <- data[data$gene == g, ]
  if (nrow(rows) == 1) {
    # Homozygous - real allele first, then placeholder
    alleles <- c(alleles, rows$allele, "Homozygous")
    normal_cov <- c(normal_cov, rows$normal, 0)
    tumor_cov <- c(tumor_cov, rows$tumor, 0)
  } else {
    alleles <- c(alleles, rows$allele)
    normal_cov <- c(normal_cov, rows$normal)
    tumor_cov <- c(tumor_cov, rows$tumor)
  }
}

mdata <- as.numeric(c(WESnormalcoverage, WEStumorcoverage, tumor_purity))

# flag low coverage alleles: if coverage of allele in normal is less than 20% of normal exome coverage or if the coverage of allele in normal is less than 5
normal_coverage_ratio <- normal_cov/WESnormalcoverage
low_normal_coverage <- ifelse(
  alleles == "Homozygous", 0,
  ifelse(normal_coverage_ratio < 0.2 | normal_cov < 5, 1, 0)
)

if (sum(low_normal_coverage, na.rm = TRUE) == sum(alleles != "Homozygous")) {
  warning("WARNING: All non-homozygous alleles are flagged as low coverage in normal. All tumor ratios will be set to 'Flag'. Check --WESnormalcoverage value and coverage file.")
}

# calculate HLA/Exome coverage ratio in tumor.
tumor_coverage_ratio <- tumor_cov/WEStumorcoverage
tumor_ratios <- ifelse(
	alleles == "Homozygous", 0,
	ifelse(low_normal_coverage == 1 , "Flag", tumor_coverage_ratio)
)
# tumor_ratios <- as.vector(tumor_ratios)

# flag if the other allele is flagged
tumor_ratios_flagged <- tumor_ratios
for(i in c(2, 4, 6)) {
  if(tumor_ratios_flagged[i-1] == "Flag" || tumor_ratios_flagged[i] == "Flag") tumor_ratios_flagged[(i-1):i] <- "Flag"
}

# limit of tumor ratio for loh
loh_tumor_ratio_limit <- ifelse(
	alleles == "Homozygous", 0,
	ifelse(tumor_ratios_flagged == "Flag" , 0, (normal_coverage_ratio * (1 - tumor_purity) + normal_coverage_ratio) / 2
))


tumor_cov_below_exp <- ifelse(
	alleles == "Homozygous", "Homozygous",
	ifelse(tumor_ratios_flagged == "Flag" , "Flag",
	ifelse(tumor_coverage_ratio < loh_tumor_ratio_limit, "Low", "Normal"
)))
# tumor_cov_below_exp <- as.data.frame(tumor_cov_below_exp)


limit_tumor_1_2 <- NULL
for(i in c(2, 4, 6)) {
  limit_tumor_1_2[i] <- normal_coverage_ratio[i] * (1 - tumor_purity) / normal_coverage_ratio[i-1]
  limit_tumor_1_2[i-1] <- normal_coverage_ratio[i-1] * (1 - tumor_purity) / normal_coverage_ratio[i]
}

# --- Warn if any limit_tumor_1_2 is NaN (division by zero from zero normal_coverage_ratio) ---
nan_limits <- which(is.nan(limit_tumor_1_2))
if (length(nan_limits) > 0) {
  warning(paste0("WARNING: limit_tumor_1_2 is NaN at position(s) ", paste(nan_limits, collapse = ", "),
                 " (allele(s): ", paste(alleles[nan_limits], collapse = ", "),
                 ") — likely due to zero normal coverage ratio. AIB calculation for these alleles will be unreliable."))
}

tumor_aib <- NULL
for(i in c(2, 4, 6)) {
  tumor_aib[i] <- ifelse(
    tumor_coverage_ratio[i] / tumor_coverage_ratio[i-1] < limit_tumor_1_2[i],
    "Low",
    "Normal"
  )
  tumor_aib[i-1] <- ifelse(
    tumor_coverage_ratio[i-1] / tumor_coverage_ratio[i] < limit_tumor_1_2[i-1],
    "Low",
    "Normal"
  )
}


aib_or_lowcov <- ifelse(
	alleles == "Homozygous", "Homozygous",
	ifelse(tumor_ratios_flagged == "Flag" , "Flag",
	ifelse(tumor_cov_below_exp == "Low" & tumor_aib == "Low", 1,
	ifelse(tumor_cov_below_exp == "Low", 0.5,
	ifelse(tumor_aib == "Low", 0.5, 0)
		))))

for(i in c(2,4,6)){
	if(alleles[i-1] == "Homozygous" || alleles[i] == "Homozygous") aib_or_lowcov[c(i-1, i)] <- "Homozygous"
}


result <- as.data.frame(cbind(alleles, aib_or_lowcov))

write.table(result, file = outname, sep = "\t", quote = F, row.names = F)
write.table(cbind(result , loh_tumor_ratio_limit, tumor_cov_below_exp), file = paste0(resultdir, "/low_cov.txt"), sep = "\t", quote = F, row.names = F)
write.table(cbind(result , limit_tumor_1_2, tumor_aib), file = paste0(resultdir, "/aib.txt"), sep = "\t", quote = F, row.names = F)

