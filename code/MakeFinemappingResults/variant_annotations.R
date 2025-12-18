#!/usr/bin/env Rscript
#' This code generates matching sets of variants for the cis-eQTL variants and the trans-eQTL variants

# Load libraries
library(tidyverse)
library(data.table)
library(arrow)

# Declare constants

# Declare function definitions


#' Match background SNPs to focal SNPs based on MAF and TSS density
#'
#' This function matches each SNP in a focal dataset to up to 1000 SNPs in a background dataset
#' based on similarity in Minor Allele Frequency (MAF) and proximity in Transcription Start Site (TSS) density.
#'
#' @param focal A data.frame or data.table containing focal SNPs. Must include columns 'MAF' and 'TSS_density'.
#' @param background A data.frame or data.table containing background SNPs. Must include columns 'MAF' and 'TSS_density'.
#'
#' @return A data.table of matched background SNPs for each focal SNP, with an added column 'focal_index'
#'         indicating the index of the focal SNP they were matched to.
#'
#' @details For each focal SNP, the function:
#' \itemize{
#'   \item Selects background SNPs with MAF within ±0.02 of the focal SNP's MAF.
#'   \item Computes the absolute difference in TSS density.
#'   \item Sorts candidates by TSS density difference.
#'   \item Returns up to 1000 closest matches.
#' }
#'
#' @examples
#' \dontrun{
#' focal_snps <- data.frame(MAF = c(0.1, 0.2), TSS_density = c(0.5, 0.7))
#' background_snps <- data.frame(MAF = runif(10000), TSS_density = runif(10000))
#' matched <- control_matching_tssd_maf(focal_snps, background_snps)
#' }
#'
#' @import data.table
#' @export
control_matching_tssd_maf_l2 <- function(focal, background, n_target=1000, n_min=0.1) {
  library(parallel)

  # Convert to data.table for speed
  focal_dt <- as.data.table(focal)
  background_dt <- as.data.table(background)

  delta_l2 <- sd(c(focal_dt[["L2"]], background_dt[["L2"]]), na.rm=T) * 0.1
  message(sprintf("Delta L2: %.2f", delta_l2))

  # Set keys for fast subsetting
  setkey(background_dt, maf)

  # Function to match one focal SNP
  match_snps <- function(focal_index, focal_l2, focal_maf, focal_tss_density, background_dt, delta_l2, delta_maf = 0.02, n = n_target) {
    candidates <- background_dt[
      focal_index != variant_index
        & maf >= (focal_maf - delta_maf)
        & maf <= (focal_maf + delta_maf)
        & L2 >= (focal_l2 - delta_l2)
        & L2 <= (focal_l2 + delta_l2)
        & tss_density == focal_tss_density]

    # Sample n with replacement
    if (nrow(candidates) >= (n_min * n_target)) {
      sampled <- candidates[sample(seq_len(nrow(candidates)), size=n, replace=T),]
    } else {
      return(data.table())
    }
    return(sampled)
  }

  cat("Starting matching...")
  cat("Starting matching...")


  match_wrapper <- function(i) {
    row <- focal_dt[i]
    cat(sprintf("\rProgress: %d/%d", i, nrow(focal_dt)))
    match_snps(row$variant_index, row$L2, row$maf, row$tss_density, background_dt, delta_l2 = delta_l2)
  }

  num_cores <- 15

  matched_list <- mclapply(1:nrow(focal_dt), match_wrapper, mc.cores = num_cores)

  names(matched_list) <- focal_dt$variant_index

  # Combine results if needed
  matched_dt <- rbindlist(matched_list, idcol = "focal_index")

  return(matched_dt)
}


calculate_tss_density <- function(variants_dt, gene_dt, radius=500000) {
  gene_dt[, chromosome := as.integer(as.character(seqid))]
  gene_dt[, .(chromosome, start = TSS, end = TSS, bp = TSS)]

  # Create a copy with start and end of the window
  variants_dt[, `:=`(start = bp - radius, end = bp + radius)]

  # Set keys for fast joins
  setkey(variants_dt, chromosome, start, end)
  setkey(gene_dt, chromosome, start, end)

  # Perform a non-equi join to count TSSs in the window
  tss_density <- foverlaps(
    x = variants_dt,
    y = gene_dt,
    by.x = c("chromosome", "start", "end"),
    by.y = c("chromosome", "start", "end"),
    type = "any", nomatch = 0L
  )[, .N, by = .(chromosome, bp)]

  # Merge back with original variants
  variants_dt <- merge(variants_dt, tss_density, by = c("chromosome", "bp"), all.x = TRUE)
  variants_dt[is.na(N), N := 0]# Replace NA with 0
  setnames(variants_dt, "N", "tss_density")

  return(variants_dt[, .(variant_index, tss_density)])
}


read_maf_table <- function(maf_path, variant_indices, sample_sizes) {
  cohort_names <- names(sample_sizes)

  # Read MAF table
  maf_table <- open_dataset(maf_path) %>% filter(variant_index %in% variant_indices) %>% collect()

  maf_matrix <- maf_table %>% as.data.table() %>% as.matrix(rownames=1)

  maf_average <- rowSums(t(t(maf_matrix[,cohort_names]) * sample_sizes[cohort_names]), na.rm=T) / rowSums(t(t(ceiling(maf_matrix[,cohort_names])) * sample_sizes[cohort_names]), na.rm=T)

  # Return MAF
  return(maf_average)
}

read_ld_scores <- function(ld_scores_path) {
  ld_scores <- bind_rows(mapply(fread, Sys.glob(ld_scores_path), SIMPLIFY=F))
  return(ld_scores)
}


# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Process input
  resources <- "resources.txt"
  variables <- deframe(read_delim(resources, delim="\t"))

  available_genes <- fread(variables["available_genes"], header=F)$V1

  available_variants <- fread(variables["available_variants"], header=F)$V1

  # Get sample sizes from master table
  sample_sizes <- fread(variables["master_table"]) %>% select(cohort_new_name, dataset_n) %>% deframe()

  hgsv_classes <- c("INS" = "Insertion", "INV" = "Inversion", "DUP" = "Duplication", "DEL" = "Deletion")

  variant_reference <- read_parquet(variables["variant_reference"]) %>%
    filter(variant_index %in% available_variants) %>%
    mutate(
      variant_class = case_when(
        str_detect(variant, "HGSV") ~ paste0(hgsv_classes[str_extract(eff_allele, "\\w{3}")], " (hgsv)"),
        str_length(eff_allele) > str_length(non_eff_allele) ~ "Insertion",
        str_length(eff_allele) < str_length(non_eff_allele) ~ "Deletion",
        str_length(eff_allele) == str_length(non_eff_allele) ~ "SNP")
    )

  gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
    filter(type == "gene") %>% select(gene_id, gene_name, seqid, start, end, strand) %>%
    mutate(TSS = if_else(strand == "+", start, end))

  # Get background set of genes
  # I.e. the maximum set of genes for which we could hypothetically detect eQTLs
  # - The intersect of the genes in the gene-reference file
  # - and in the set of ensembl IDs for which the meta-analysis was run
  gene_reference_filtered <- gene_reference %>% filter(gene_id %in% available_genes)

  tss_density_table <- calculate_tss_density(variants_dt = as.data.table(variant_reference), as.data.table(gene_reference_filtered), radius=500000)
  maf_table <- read_maf_table(maf_path = variables["maf_file"], variant_reference$variant_index, sample_sizes)

  maf_table <- enframe(maf_table, name = "variant_index", value = "maf") %>%
    mutate(variant_index = as.integer(variant_index))

  ld_scores <- read_ld_scores(variables["ld_scores_1kg"]) %>% select(SNP, L2)

  variant_annotations <- left_join(variant_reference, tss_density_table, by = "variant_index") %>%
    left_join(maf_table, by = "variant_index") %>%
    left_join(ld_scores, by = c("variant" = "SNP"))

  variant_annotations.parquet <- "variant_annotations.parquet"
  write_parquet(variant_annotations, variant_annotations.parquet)
  rm(maf_table, ld_scores, tss_density_table)
  variant_annotations <- read_parquet(variant_annotations.parquet)
  # Process output
  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type))

  annotated_finemapped_variants <- left_join(
    variant_annotations,
    independent_variant_set_annot_filtered %>%
      distinct(variant_index, typeCisTrans))

  # Generate matching set for just the cis-eQTLs and for just the trans-eQTLs
  matching_set_cis <- control_matching_tssd_maf_l2(
    focal = annotated_finemapped_variants %>% filter(typeCisTrans == "cis", !is.na(L2)),
    background = annotated_finemapped_variants %>% filter(is.na(typeCisTrans)), n_target=100)

  fwrite(matching_set_cis, "matching_set_cis_eQTLs_20251017.txt.gz", sep="\t", row.names=F, col.names=T)
  write_parquet(matching_set_cis, "matching_set_cis_eQTLs_20251017.parquet")

  matching_set_trans <- control_matching_tssd_maf_l2(
    focal = annotated_finemapped_variants %>% filter(typeCisTrans == "trans", !is.na(L2)),
    background = annotated_finemapped_variants %>% filter(is.na(typeCisTrans)))

  fwrite(matching_set_trans, "matching_set_trans_eQTLs_20251017.txt.gz", sep="\t", row.names=F, col.names=T)
  write_parquet(matching_set_trans, "matching_set_trans_eQTLs_20251017.parquet")

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}