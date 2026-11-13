 result <- lapply(1:nrow(hla_la_all), function(z){

 la <- unlist(hla_la_all[z,2:7])
 hd <- unlist(hla_hd_all[z,2:7])

# named vectors of hla alleles
la <- setNames(la, paste0("V", 1:6))
hd <- setNames(hd, paste0("V", 1:6))

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

 return(genomic_alleles_new)})



 result1 <- lapply(1:length(result), function(x) c(unlist(result[[x]][1,2:3], use.names=F), unlist(result[[x]][2,2:3], use.names=F), unlist(result[[x]][3,2:3],use.names=F) ) )
result2 <- do.call('rbind', result1)
rownames(result2) <- rownames(hla_la_all)
 write.xlsx(result2, file='test_out_04302026__.xlsx')