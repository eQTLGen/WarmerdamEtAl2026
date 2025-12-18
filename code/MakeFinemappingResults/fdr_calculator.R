#!/usr/bin/env Rscript
#' Code to calculate FDR for cis- and trans- eQTLs

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

  available_genes <- fread(variables["available_genes"], header = FALSE)$V1

  gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
    filter(type == "gene") %>%
    select(gene_id, gene_name, seqid, start, end, strand) %>%
    mutate(TSS = if_else(strand == "+", start, end))

  gene_reference_filtered <- gene_reference %>%
    filter(gene_id %in% available_genes)

  tested_genes <- gene_reference %>%
    filter(gene_id %in% available_genes) %>%
    pull(gene_id)

  available_genes_cis <- gene_reference %>%
    filter(
      seqid != "X",
      gene_id %in% available_genes,
      !(seqid == 6 & start <= 34000000 & 25000000 <= end)
    ) %>%
    pull(gene_id)

  meta_analysed <- fread(variables["independent_variant_set"]) %>%
    select(phenotype, variant_index, beta, standard_error, type, p_value, pip, locus_str)

  nr_of_primary_cis_eqtl_genes <- meta_analysed %>%
    filter(phenotype %in% available_genes_cis) %>%
    group_by(phenotype, locus_str) %>%
    slice_max(abs(beta/standard_error), with_ties=F) %>%
    ungroup() %>%
    filter(type == "cis") %>%
    distinct(phenotype) %>%
    nrow()

  nr_of_primary_trans_eqtls <- meta_analysed %>%
    filter(phenotype %in% tested_genes) %>%
    group_by(phenotype, locus_str) %>%
    slice_max(abs(beta/standard_error), with_ties=F) %>%
    ungroup() %>%
    filter(type != "cis") %>%
    distinct(phenotype, locus_str) %>%
    nrow()

  genome_indep_variants <- 2158129
  genome_size <- 3.1e9
  window_size <- 2e6
  num_genes <- length(tested_genes)
  num_genes_cis <- length(available_genes_cis)
  p_threshold <- 5e-8
  observed_significant <- nr_of_primary_cis_eqtl_genes

  fraction_per_gene <- window_size / genome_size
  indep_variants_per_gene <- genome_indep_variants * fraction_per_gene
  expected_fp_cis <- num_genes_cis * indep_variants_per_gene * p_threshold
  fdr_cis <- expected_fp_cis / observed_significant * 100

  cat("=== CIS-eQTL ===\n")
  cat("Nr of tested genes:", num_genes_cis, "\n")
  cat("Independent variants per gene:", round(indep_variants_per_gene), "\n")
  cat("Expected false positives:", round(expected_fp_cis, 2), "\n")
  cat("FDR (%):", round(fdr_cis, 4), "\n\n")

  observed_trans <- nr_of_primary_trans_eqtls
  num_scans <- num_genes
  expected_fp_trans <- num_scans * genome_indep_variants * p_threshold
  fdr_trans <- expected_fp_trans / observed_trans * 100

  cat("=== TRANS-eQTL ===\n")
  cat("Nr of primary trans-eQTLs:", observed_trans, "\n")
  cat("Expected false positives:", round(expected_fp_trans), "\n")
  cat("FDR (%):", round(fdr_trans, 2), "\n")

  # Work
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}