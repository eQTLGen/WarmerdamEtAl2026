#!/usr/bin/env Rscript
#` This code summarises output of the gwas gene prioritization

# Load libraries
library(tidyverse)
library(data.table)
library(readxl)

# Declare constants

# Declare function definitions

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))

  gwas_selection <- readxl::read_xlsx(variables["eqtlgen_gwass"]) %>%
    select(trait_identifier, trait_name)
  # Fill in!
  output_prefix <- variables["network_analyses_prefix"]

  # Process input
  gene_prioritizations_all <- fread(sprintf("%s_prioritization_perm_results.rds", output_prefix)) %>%
    mutate(gene_score_p_adj_b = p.adjust(gene_score_p, method = "bonferroni"), gene_score_p_adj_bh = p.adjust(gene_score_p, method = "BH"),
           fdr_significant = gene_score_p_adj_bh < 0.05, Bonferroni_significant = gene_score_p_adj_b < 0.05) %>%
    select(gwas_id, gene, symbol, gene_score, gene_score_mean, gene_score_sd, gene_score_z, gene_score_p, gene_score_p_adj_b, gene_score_p_adj_bh) %>%
    left_join(gwas_selection, by=c("gwas_id" = "trait_identifier"))

  fwrite(gene_prioritizations_all, "gene_prioritization_results_all.txt", sep="\t", row.names=F, col.names=T)

  significant_results <- gene_prioritizations_all %>%
    filter(gene_score_p_adj_bh < 0.05, abs(gene_score) > 1)

  significant_results %>%
    group_by(sign(gene_score)) %>%
    summarise(n = n()) %>%
    mutate(p = n / sum(n) * 100)

  significant_results %>%
    group_by(gwas_id) %>%
    summarise(n=n()) %>%
    summarise(n_studies = n(), sum = sum(n), mean = sum / n_studies, mean_all = sum / 70)
  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}