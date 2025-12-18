#!/usr/bin/env Rscript
#' Script that extracts finemapped eQTLs from the INTERVAL sQTL data

# Load libraries
library(tidyverse)
library(data.table)
library(arrow)


# Declare constants


# Define functions


# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Process input

  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))
  harmonized_sQTL_dataset <- "" # The harmonized dataset with sQTLs
  phenotype_summary_path <- "" # The phenotype summary from INTERVAL

  variables <- deframe(read_delim(resources, delim = "\t"))
  location <- str_split(variables["available_genes"], "freeze3")[[1]][1]

  available_genes <- fread(variables["available_genes"], header = FALSE)$V1

  sqtl_ds <- open_dataset(harmonized_sQTL_dataset) # The harmonized dataset with sQTLs

  independent_set_finemapped_cis <- fread(variables["independent_variant_set"]) %>%
    filter(type == "cis") %>% select(variant_index, phenotype, chromosome)

  phenotype_summary <- fread(phenotype_summary_path) %>% #
    distinct(gene_id, phenotype_id) %>%
    separate_longer_delim(gene_id, delim=",")

  independent_set_finemapped_cis_with_phenotype_ids_per_chr <- independent_set_finemapped_cis %>%
    left_join(phenotype_summary, by = c("phenotype" = "gene_id"), relationship = "many-to-many") %>%
    group_by(chromosome) %>%
    group_split(.keep = TRUE)

  queried_sqtls <- mapply(
    function(potential_cis_sqtls) {
      chromosome <- potential_cis_sqtls$chromosome[1]
      cat("Chromosome:", chromosome)
      sqtl_chrom <- sqtl_ds %>%
        filter(CHR == !!chromosome) %>% collect()

      sqtl_chrom %>%
        right_join(potential_cis_sqtls, by = c("phenotype_id", "variant_index")) %>%
        collect()
    }
  , independent_set_finemapped_cis_with_phenotype_ids_per_chr, SIMPLIFY=F, USE.NAMES = T)

  queried_sqtls_bound <- bind_rows(queried_sqtls)

  sqtl_stats_per_cis_eqtl <- queried_sqtls_bound %>%
    mutate(sqtl_z = beta/se) %>%
    group_by(variant_index, phenotype) %>%
    summarise(sqtl_phenotypes = paste0(phenotype_id, collapse=","),
              sqtl_top_phenotype = ifelse(any(!is.na(sqtl_z)), phenotype_id[which.max(abs(sqtl_z))], NA),
              sqtl_z_scores = paste0(sqtl_z, collapse=","),
              sqtl_max_abs_z_score = ifelse(any(!is.na(sqtl_z)), max(abs(sqtl_z), na.rm=T), NA),
              sqtl_mean_abs_z_score = ifelse(any(!is.na(sqtl_z)), mean(abs(sqtl_z), na.rm=T), NA),
              sqtl_mean_abs_z_score = ifelse(any(!is.na(sqtl_z)), mean(abs(sqtl_z), na.rm=T), NA),
    )


  independent_set_finemapped_cis <- fread(variables["independent_variant_set"]) %>%
    filter(type == "cis") %>% select(variant_index, variant, phenotype, chromosome) %>%
    left_join(sqtl_stats_per_cis_eqtl, by = c("variant_index", "phenotype"))


  fwrite(independent_set_finemapped_cis, "sqtls_lookup_for_all_cis_eqtls.txt.gz", sep="\t", row.names=F, col.names=T, quote=F)

  # Work
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}