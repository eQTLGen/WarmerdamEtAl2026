#!/usr/bin/env Rscript


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
  # Process input

  # Define variables
  outcome_variables <- list()

  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))

  available_genes <- fread(variables["available_genes"], header=F)$V1
  all_significant <- fread(variables["all_significant"])

  gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
    filter(type == "gene") %>% select(gene_id, gene_name, seqid, start, end, strand) %>%
    mutate(TSS = if_else(strand == "+", start, end))

  finemapping_log <- fread(variables["finemapping_log"]) %>%
    select(phenotype = gene, chromosome, locus_start = start, locus_end = end) %>% as.data.table()

  # Load the main file of independent variants
  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type))

  # Work
  loci <- independent_variant_set_annot_filtered %>%
    distinct(phenotype, chromosome, locus_start, locus_end, locus_str) %>% as.data.table()

  setkey(loci, phenotype, chromosome, locus_start, locus_end)
  setkey(finemapping_log, phenotype, chromosome, locus_start, locus_end)

  finemapped_loci <- foverlaps(
    loci,
    finemapping_log,
    by.x = c("phenotype", "chromosome", "locus_start", "locus_end"),
    by.y = c("phenotype", "chromosome", "locus_start", "locus_end"),
    type = "any",
    nomatch = 0
  )

  fwrite(finemapped_loci, "finemapped_eqtls_locus_str_to_original_finemapping_boundaries_per_cluster.txt.gz", sep="\t", row.names=F, col.names=T)

  finemapping_loci_summarised <- finemapped_loci %>%
    group_by(phenotype, chromosome, locus_str) %>%
    summarise(locus_start = min(locus_start), locus_end = max(locus_end))

  fwrite(finemapping_loci_summarised, "finemapped_eqtls_locus_str_to_original_finemapping_boundaries_summarised.txt.gz", sep="\t", row.names=F, col.names=T)

  finemapping_loci_summarised_dt <- finemapping_loci_summarised %>% select(phenotype, chromosome, locus_start, locus_end, locus_str) %>% as.data.table()
  setkey(finemapping_loci_summarised_dt, phenotype, chromosome, locus_start, locus_end)

  all_significant[, `:=`(locus_start = bp, locus_end = bp)]
  setkey(all_significant, phenotype, chromosome, locus_start, locus_end)

  result <- foverlaps(
    all_significant,
    finemapping_loci_summarised_dt,
    by.x = c("phenotype", "chromosome", "locus_start", "locus_end"),
    by.y = c("phenotype", "chromosome", "locus_start", "locus_end"),
    type = "within",
    nomatch = 0
  )

  cis_window <- 1e6

  summarised <- result %>%
    group_by(phenotype, chromosome, locus_start, locus_end, locus_str) %>%
    summarise(n_variants = n(), min_p_val = min(p_value), cis = cis_1e6[which.min(p_value)], bp = bp[which.min(p_value)], variant_index=variant_index[which.min(p_value)])

  summarised %>% filter(min_p_val < (5e-8)) %>%
    group_by(cis) %>%
    tally()

  summarised %>% filter(min_p_val < (5e-8 / 16742)) %>%
    group_by(cis) %>%
    tally()
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}