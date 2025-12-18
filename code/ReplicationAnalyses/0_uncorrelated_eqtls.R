#!/usr/bin/env Rscript
#' This code makes a list of eQTLs that is uncorrelated.
#' The output is used for calculating Rbs

# Load libraries
library(data.table)
library(arrow)
library(tidyverse)
library(dtplyr)
library(dbplyr)



prune_susie_output <- function(chromosome_eqtl_table, ld_func, r2_threshold) {
  chromosome_eqtl_table <- chromosome_eqtl_table %>%
    group_by(variant_index) %>%
    slice_max(max_lbf)

  if (nrow(chromosome_eqtl_table) == 1) {
    return(chromosome_eqtl_table)
  }

  # Get the LD matrix
  ld_matrix <- ld_func(chromosome_eqtl_table$chromosome[1], chromosome_eqtl_table$variant_index)
  ld_matrix <- ld_matrix[
    as.character(chromosome_eqtl_table$variant_index),
    as.character(chromosome_eqtl_table$variant_index)]

  if (!all(rownames(ld_matrix) == as.character(chromosome_eqtl_table$variant_index))) {
    stop("LD matrix row names must match variant IDs in variants table")
  }

  # Create new columns for tagging information
  chromosome_eqtl_table$tagging_variants <- NA
  chromosome_eqtl_table$tagging_r <- NA

  # Sort variants by max_lbf (descending)
  chromosome_eqtl_table <- chromosome_eqtl_table %>% ungroup() %>% arrange(desc(max_lbf), desc(z_score_abs))

  # Keep track of removed variants
  to_remove <- rep(FALSE, nrow(chromosome_eqtl_table))
  names(to_remove) <- rownames(ld_matrix)

  # Iterate through variants
  for (index_variant in as.character(chromosome_eqtl_table$variant_index)) {
    if (to_remove[index_variant]) next  # Skip already tagged variants

    ld_values <- ld_matrix[index_variant, ]  # Get r^2 values

    # Identify variants in high LD
    which_variants_high_ld <- which(
      (ld_values^2) > r2_threshold &
        names(ld_values) != index_variant)
    high_ld_variants <- ld_values[which_variants_high_ld]

    if (length(high_ld_variants) > 0) {
      # Tag these variants
      chromosome_eqtl_table$tagging_variants[chromosome_eqtl_table$variant_index == as.integer(index_variant)] <- paste0(names(high_ld_variants), collapse=",")
      chromosome_eqtl_table$tagging_r[chromosome_eqtl_table$variant_index == as.integer(index_variant)] <- paste0(unname(high_ld_variants), collapse=",")

      # Mark them for removal
      to_remove[names(high_ld_variants)] <- TRUE
    }
  }

  return(chromosome_eqtl_table)
}



# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Define variables
  resources <- "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/users/umcg-rwarmerdam/resources.txt"
  variables <- deframe(read_delim(resources, delim="\t"))
  location <- str_split(variables["available_genes"], "freeze3")[[1]][1]

  # Load the main file of independent variants
  eqtl_table_raw <- fread(variables["independent_variant_set"])


  eqtl_table <- eqtl_table_raw %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type)) %>%
    select(phenotype, chromosome, variant, variant_index, bp, max_lbf, pip, beta, standard_error, typeCisTrans, type) %>%
    mutate(z_score_abs = abs(beta / standard_error))

  precalculated_ld_path <- "../ld_matrices_20250509.rds"
  if (file.exists(precalculated_ld_path)) {
    ld_matrices <- readRDS(precalculated_ld_path)
  }

  ld_func_mat <- function(chromosome, variant_indices) {
    chromosome_ld_matrix <- ld_matrices[[as.character(chromosome)]]
    return(chromosome_ld_matrix[as.character(variant_indices), as.character(variant_indices)])
  }

  output_prune_function <- bind_rows(mapply(function(.x) {
    prune_susie_output(.x, ld_func_mat, 0.1)
  }, eqtl_table %>%
                               group_by(chromosome, phenotype) %>%
                               group_split(), SIMPLIFY=F))

  output_pruning_with_flag <- output_prune_function %>%
    mutate(across(c(tagging_variants, tagging_r)
      , ~ lapply(str_split(.x, ","), as.numeric))) %>%
    group_by(chromosome, phenotype) %>%
    mutate(high_r2_flag = variant_index %in% c(unlist(tagging_variants))) %>%
    select(variant_index, phenotype, stringent_r2_flag = high_r2_flag)

  eqtl_table_raw_with_01 <- eqtl_table_raw %>% inner_join(output_pruning_with_flag) %>%
    filter(!stringent_r2_flag)
  fwrite(eqtl_table_raw_with_01, "independent_eQTLs_stringently_clumped_r2_01_for_Rb_20251106.txt.gz", sep="\t", row.names=F, col.names=T)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}