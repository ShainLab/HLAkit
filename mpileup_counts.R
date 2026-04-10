#!/usr/bin/env Rscript
#!/opt/homebrew/bin Rscript
#!/opt/homebrew/share/man/man1 Rscript.1

# Map back the normal and tumor mpileup output to the HLA somatic mutations file
library(optparse)

option_list <- list( 
    make_option(c("-s", "--somatic_mutations"), 
        type = "character",
        help = "somatic mutations file"
        ),
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
somatic_mutations <- opt$somatic_mutations
normal_mpileupout <- opt$normal_mpileupout
tumor_mpileupout_mapq0 <- opt$tumor_mpileupout_mapq0
tumor_mpileupout_mapqnonzero <- opt$tumor_mpileupout_mapqnonzero

# --- Argument validation ---
if (is.null(somatic_mutations)) stop("ERROR: -s/--somatic_mutations is required but not provided.")
if (is.null(normal_mpileupout)) stop("ERROR: -n/--normal_mpileupout is required but not provided.")
if (is.null(tumor_mpileupout_mapq0)) stop("ERROR: -m/--tumor_mpileupout_mapq0 is required but not provided.")
if (is.null(tumor_mpileupout_mapqnonzero)) stop("ERROR: -u/--tumor_mpileupout_mapqnonzero is required but not provided.")

check_file <- function(f) {
    if (!file.exists(f)) stop(paste("File not found:", f))
    if (file.info(f)$size == 0) stop(paste("File is empty:", f))
}

check_file(somatic_mutations)
check_file(normal_mpileupout)
check_file(tumor_mpileupout_mapq0)
check_file(tumor_mpileupout_mapqnonzero)

cat("Input files check complete.\n")

#output filename with mpileup counts for both normal and tumor
outfile <- gsub(".txt", ".mpileup_counts.txt", somatic_mutations)

# read files
somatic_mutations <- read.delim(somatic_mutations, header = T, sep = '\t', quote = "", stringsAsFactors = F, colClasses = "character")
normal_mpileupout <- read.delim(normal_mpileupout, header = F, sep = '\t', quote = "", stringsAsFactors = F)
tumor_mpileupout_mapq0 <- read.delim(tumor_mpileupout_mapq0, header = F, sep = '\t', quote = "", stringsAsFactors = F)
tumor_mpileupout_mapqnonzero <- read.delim(tumor_mpileupout_mapqnonzero, header = F, sep = '\t', quote = "", stringsAsFactors = F)

if (nrow(somatic_mutations) == 0) stop("ERROR: Somatic mutations file has no data rows.")

required_cols <- c("REF", "ALT")
missing_cols <- setdiff(required_cols, colnames(somatic_mutations))
if (length(missing_cols) > 0) stop(paste("ERROR: Somatic mutations file is missing required column(s):", paste(missing_cols, collapse = ", ")))

if (nrow(normal_mpileupout) == 0) stop("ERROR: Normal mpileup output file has no data rows.")
if (nrow(tumor_mpileupout_mapq0) == 0) stop("ERROR: Tumor MAPQ0 mpileup output file has no data rows.")
if (nrow(tumor_mpileupout_mapqnonzero) == 0) stop("ERROR: Tumor MAPQ-nonzero mpileup output file has no data rows.")

#all mpileup files must match somatic mutations
n_mut <- nrow(somatic_mutations)
if (nrow(normal_mpileupout) != n_mut) {
  stop(paste0("ERROR: Row count mismatch — somatic mutations has ", n_mut, " rows but normal mpileup has ",
              nrow(normal_mpileupout), " rows. Mpileup files must have one row per mutation."))
}
if (nrow(tumor_mpileupout_mapq0) != n_mut) {
  stop(paste0("ERROR: Row count mismatch — somatic mutations has ", n_mut, " rows but tumor MAPQ0 mpileup has ",
              nrow(tumor_mpileupout_mapq0), " rows. Mpileup files must have one row per mutation."))
}
if (nrow(tumor_mpileupout_mapqnonzero) != n_mut) {
  stop(paste0("ERROR: Row count mismatch — somatic mutations has ", n_mut, " rows but tumor MAPQ-nonzero mpileup has ",
              nrow(tumor_mpileupout_mapqnonzero), " rows. Mpileup files must have one row per mutation."))
}

# check mpileup files columns
for (nm in c("normal_mpileupout", "tumor_mpileupout_mapq0", "tumor_mpileupout_mapqnonzero")) {
  df <- get(nm)
  if (ncol(df) < 5) stop(paste0("ERROR: ", nm, " has only ", ncol(df), " column(s); expected at least 5 (chrom, pos, ref, depth, bases). Check mpileup output format."))
}

# index <- match(mpileup_pos, maf_pos)
refbase <- somatic_mutations$REF
altbase <- somatic_mutations$ALT

result <- NULL

#find indels
indels <- which((nchar(refbase) != nchar(altbase)) | refbase == "-" | refbase == "." | altbase == "-" | altbase == "." )
snvs <- which(nchar(refbase) == nchar(altbase))

message(paste0("NOTE: ", length(indels), " indel(s) and ", length(snvs), " SNV/DNP(s) identified."))

# find ref and alt counts for indels
if(length(indels) != 0){
    indel_normal_mpileupout <- normal_mpileupout[indels, ]
    indel_tumor_mpileupout_mapq0 <- tumor_mpileupout_mapq0[indels, ]
    indel_tumor_mpileupout_mapqnonzero <- tumor_mpileupout_mapqnonzero[indels, ]
    indel_somatic_mutations <- somatic_mutations[indels, ]

    snv_normal_mpileupout <- normal_mpileupout[-(c(indels)), ]
    snv_tumor_mpileupout_mapq0 <- tumor_mpileupout_mapq0[-(c(indels)), ]
    snv_tumor_mpileupout_mapqnonzero <- tumor_mpileupout_mapqnonzero[-(c(indels)), ]
    snv_somatic_mutations <- somatic_mutations[-(c(indels)), ]

    normal_indel_counts <- lapply(1:nrow(indel_normal_mpileupout), function(x){
        if(indel_normal_mpileupout[x,4] > 0){
            read <- strsplit(indel_normal_mpileupout[x, 5], split = "")[[1]]
            ref <- sum(length(grep("\\.", read)), length(grep(",", read)))
            alt <- sum(length(grep("\\+", read)), length(grep("\\*", read)))
        }
        else{
            ref <- 0
            alt <- 0
        }
        return(c(ref, alt))
        })
    normal_indel_counts <- do.call('rbind', normal_indel_counts)

    tumor_indel_counts_mapq0 <- lapply(1:nrow(indel_tumor_mpileupout_mapq0), function(x){
        if(indel_tumor_mpileupout_mapq0[x,4] > 0){
            read <- strsplit(indel_tumor_mpileupout_mapq0[x, 5], split = "")[[1]]
            ref <- sum(length(grep("\\.", read)), length(grep(",", read)))
            alt <- sum(length(grep("\\+", read)), length(grep("\\*", read)))
        }
        else{
            ref <- 0
            alt <- 0
        }
        return(c(ref, alt))
        })
    tumor_indel_counts_mapq0 <- do.call('rbind', tumor_indel_counts_mapq0)

    tumor_indel_counts_mapqnonzero <- lapply(1:nrow(indel_tumor_mpileupout_mapqnonzero), function(x){
        if(indel_tumor_mpileupout_mapqnonzero[x,4] > 0){
            read <- strsplit(indel_tumor_mpileupout_mapqnonzero[x, 5], split = "")[[1]]
            ref <- sum(length(grep("\\.", read)), length(grep(",", read)))
            alt <- sum(length(grep("\\+", read)), length(grep("\\*", read)))
        }
        else{
            ref <- 0
            alt <- 0
        }
        return(c(ref, alt))
        })
    tumor_indel_counts_mapqnonzero <- do.call('rbind', tumor_indel_counts_mapqnonzero)


    #write results
    indel_result <- cbind(indel_somatic_mutations, 'Normal_Ref' = 0, 'Normal_Mut' = 0, "Tumor_Ref" = 0, "Tumor_Mut" = 0)
    indel_result[, 'Tumor_Ref'] <- tumor_indel_counts_mapq0[,1]/2 + tumor_indel_counts_mapqnonzero[,1]
    indel_result[, 'Tumor_Mut'] <- tumor_indel_counts_mapq0[,2] + tumor_indel_counts_mapqnonzero[,2]
    indel_result[, c('Normal_Ref', 'Normal_Mut')] <- normal_indel_counts
    result <- indel_result

} else{
        snv_normal_mpileupout <- normal_mpileupout
        snv_tumor_mpileupout_mapq0 <- tumor_mpileupout_mapq0
        snv_tumor_mpileupout_mapqnonzero <- tumor_mpileupout_mapqnonzero
        snv_somatic_mutations <- somatic_mutations
        refbase <- snv_somatic_mutations$REF
        altbase <- snv_somatic_mutations$ALT
}

if(length(snvs) > 0){
    # find ref and alt counts for snvs
    normal_snv_counts<- lapply(1:nrow(snv_normal_mpileupout), function(x){
        if(snv_normal_mpileupout[x,4] > 0){
        read <- strsplit(snv_normal_mpileupout[x, 5], split = "")[[1]]
        ref <- sum(length(grep("\\.", read)), length(grep(",", read)))
        alt <- sum(length(grep(substring(altbase[x],1,1), read, ignore.case = T)))
        }
        else{
            ref <- 0
            alt <- 0
        }
        return(c(ref, alt))
        })
    normal_snv_counts <- do.call('rbind', normal_snv_counts)

    tumor_snv_counts_mapq0 <- lapply(1:nrow(snv_tumor_mpileupout_mapq0), function(x){
        if(snv_tumor_mpileupout_mapq0[x,4] > 0){
        read <- strsplit(snv_tumor_mpileupout_mapq0[x, 5], split = "")[[1]]
        ref <- sum(length(grep("\\.", read)), length(grep(",", read)))
        alt <- sum(length(grep(substring(altbase[x],1,1), read, ignore.case = T)))
        }
        else{
            ref <- 0
            alt <- 0
        }
        return(c(ref, alt))
        })
    tumor_snv_counts_mapq0 <- do.call('rbind', tumor_snv_counts_mapq0)

    tumor_snv_counts_mapqnonzero <- lapply(1:nrow(snv_tumor_mpileupout_mapqnonzero), function(x){
        if(snv_tumor_mpileupout_mapqnonzero[x,4] > 0){
        read <- strsplit(snv_tumor_mpileupout_mapqnonzero[x, 5], split = "")[[1]]
        ref <- sum(length(grep("\\.", read)), length(grep(",", read)))
        alt <- sum(length(grep(substring(altbase[x],1,1), read, ignore.case = T)))
        }
        else{
            ref <- 0
            alt <- 0
        }
        return(c(ref, alt))
        })
    tumor_snv_counts_mapqnonzero <- do.call('rbind', tumor_snv_counts_mapqnonzero)

    #write results
    snv_result <- cbind(snv_somatic_mutations, 'Normal_Ref' = 0, 'Normal_Mut' = 0, "Tumor_Ref" = 0, "Tumor_Mut" = 0)
    snv_result[, 'Tumor_Ref'] <- tumor_snv_counts_mapq0[,1]/2 + tumor_snv_counts_mapqnonzero[,1]
    snv_result[, 'Tumor_Mut'] <- tumor_snv_counts_mapq0[,2] + tumor_snv_counts_mapqnonzero[,2]
    snv_result[, c('Normal_Ref', 'Normal_Mut')] <- normal_snv_counts

    result <- rbind(result, snv_result)
}

if (is.null(result) || nrow(result) == 0) stop("ERROR: Result is empty after processing both indels and SNVs. Check that REF/ALT columns contain valid values and that mpileup row counts match the somatic mutations file.")

result$Tumor_MAF <- as.numeric(result$Tumor_Mut) / (as.numeric(result$Tumor_Ref) + as.numeric(result$Tumor_Mut))
result$Tumor_MAF[which(is.nan(result$Tumor_MAF))] <- 0

n_zero_depth <- sum(as.numeric(result$Tumor_Ref) + as.numeric(result$Tumor_Mut) == 0)
if (n_zero_depth > 0) warning(paste0("WARNING: ", n_zero_depth, " mutation(s) have zero tumor depth (Tumor_Ref + Tumor_Mut == 0); Tumor_MAF set to 0 for these. Check mpileup coverage at these positions."))

# result$Artifacts <- ""
# result$Artifacts <- ifelse(result$Normal_Mut > 2, "artifact",
#     ifelse(result$Tumor_Mut == 0, "artifact",
#         ifelse(result$Tumor_Mut + result$Tumor_Ref == 0, "artifact", "")))

# n_artifacts <- sum(result$Artifacts == "artifact")
# if (n_artifacts > 0) message(paste0("NOTE: ", n_artifacts, " mutation(s) flagged as artifacts (Normal_Mut > 2, or zero tumor reads)."))
# if (n_artifacts == nrow(result)) warning("WARNING: All mutations flagged as artifacts. Check mpileup input files and mutation coordinates.")


write.table(result, file = outfile, sep = "\t", quote = F, row.names = F, col.names = T)
cat("Counting Ref and Mut reads done!\n")



cat("Counting Ref and Mut reads done!\n")

