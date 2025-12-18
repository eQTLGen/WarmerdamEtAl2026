#!/usr/bin/env Rscript
#' Count the number of cis-eQTLs with a trans-eQTL

# Load libraries
library(tidyverse)
library(data.table)


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

  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))

  # Process input
  meta_analysed <- fread(variables["independent_variant_set"]) %>%
    select(phenotype, variant_index, beta, standard_error, type, p_value, pip)

  nr_of_cis_effects_potential <- meta_analysed %>%
    filter(type == "cis", !is.na(pip)) %>%
    nrow()

  nr_of_trans_effects_potential <- meta_analysed %>%
    filter(type == "trans", !is.na(pip)) %>%
    nrow()

  # Process input
  lbf_based_cis_trans_coloc <- fread(variables["lbf_based_cis_trans_coloc"])

  nrow(lbf_based_cis_trans_coloc)

  nr_of_colocalizing_cis_effects <- lbf_based_cis_trans_coloc %>%
    distinct(cis_gene, signal1) %>%
    nrow()

  nr_of_colocalizing_cis_effects
  nr_of_cis_effects_potential
  nr_of_colocalizing_cis_effects / nr_of_cis_effects_potential

  nr_of_colocalizing_trans_effects <- lbf_based_cis_trans_coloc %>%
    distinct(trans_gene, signal2) %>%
    nrow()

  nr_of_colocalizing_trans_effects
  nr_of_trans_effects_potential
  nr_of_colocalizing_trans_effects / nr_of_trans_effects_potential * 100
  # Work
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}