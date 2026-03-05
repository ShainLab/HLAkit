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

# index <- match(mpileup_pos, maf_pos)
refbase <- somatic_mutations$REF
altbase <- somatic_mutations$ALT

result <- NULL

#find indels
indels <- which((nchar(refbase) != nchar(altbase)) | refbase == "-" | refbase == "." | altbase == "-" | altbase == "." )
snvs <- which(nchar(refbase) == nchar(altbase))
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

if(length(snvs) != 0){
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


# artifacts
# if normal has >2 counts OR if normal has zero coverage OR if tumor has <4 counts OR if tumor has zero coverage, call artifact
artifacts <- unique(c(
    which(result[,"Normal_Ref"] + result[,"Normal_Mut"] == 0),
    which(result[,"Normal_Mut"] > 2),
    which(result[,"Tumor_Mut"] < 4),
    which(result[,"Tumor_Ref"] + result[,"Tumor_Mut"] == 0)))

result$Tumor_MAF <- as.numeric(result$Tumor_Mut) / (as.numeric(result$Tumor_Ref) + as.numeric(result$Tumor_Mut))
result$Artifacts <- ""

if(length(artifacts) != 0) result[artifacts, 'Artifacts'] <- 'artifact'
write.table(result, file = outfile, sep = "\t", quote = F, row.names = F, col.names = T)

cat("Counting Ref and Mut reads done!\n")

