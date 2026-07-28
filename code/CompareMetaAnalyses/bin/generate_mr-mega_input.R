#!/usr/bin/env Rscript
#' Generate plots for replication results

# Load libraries
library(tidyverse)
library(arrow)
library(data.table)
library(ggrepel)
library(patchwork)
library(argparse)

# Declare constants
old <- theme_set(theme_classic(base_size = 8))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27 / 96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  panel.grid.major = element_line(colour="grey90"),
  strip.background = element_rect(colour = NA, fill = NA),
  text = element_text(family = "Helvetica"),
  axis.title = element_text(colour = "black"),
  title = element_text(colour = "#595A5C", face = "plain")
)

# -------------------------
# Argparse CLI definition
# -------------------------
parser <- ArgumentParser(
  description = "Generate MR-MEGA input summary statistics per cohort"
)

parser$add_argument(
  "--resources", "-r",
  help = "Path to resources.txt (tab-delimited key/value table)",
  default = "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/users/umcg-rwarmerdam/resources.txt"
)

parser$add_argument(
  "--per_cohort_extractions", "-p",
  help = "Directory containing per-cohort extraction results (Arrow dataset root)",
  default = "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/ExtractPerCohortEstimates/results"
)


# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  args <- parser$parse_args(argv)

  resources <- args$resources
  per_cohort_extractions <- args$per_cohort_extractions

  variables <- deframe(read_delim(resources, delim="\t"))
  location <- str_split(variables["available_genes"], "freeze3")[[1]][1]

  #cohort_info <- fread("cohort_info.txt")
  master_table <- fread(variables["master_table"]) %>% select(cohort, cohort_new_name, dataset_n, dataset_platform, dataset_ancestry)
  sample_sizes <- master_table %>% select(cohort, dataset_n) %>% deframe()


  # Process input
  variant_reference <- read_parquet(variables["variant_reference"])

  per_cohort_summary_statistics <- open_dataset(per_cohort_extractions) %>% collect()
  per_cohort_summary_statistics <- per_cohort_summary_statistics %>%
    inner_join(variant_reference, by = "variant_index")

  variants_to_filter <- unique(per_cohort_summary_statistics$variant_index)
  all_cohort_names <- unique(per_cohort_summary_statistics$cohort)

  maf_table <- read_parquet(variables["maf_file"]) %>%
    filter(variant_index %in% variants_to_filter) %>% collect()
  maf_table_long <- maf_table %>%
    pivot_longer(cols=-variant_index, names_to="cohort_new_name", values_to = "eaf") %>%
    inner_join(master_table, by = c("cohort_new_name")) %>%
    select(variant_index, cohort, eaf)

  all_cohort_names %in% colnames(maf_table)
  all(all_cohort_names %in% maf_table_long$cohort)
  all(unique(maf_table_long$cohort) %in% all_cohort_names)

  per_cohort_summary_statistics_with_maf <- per_cohort_summary_statistics %>%
    inner_join(maf_table_long, by = c("cohort", "variant_index")) %>%
    mutate(sample_size = sample_sizes[cohort])

  # Need to prepare these columns:

  # 1) MARKERNAME – snp name
  #
  # 2) EA – effect allele
  #
  # 3) NEA – non effect allele
  #
  # 4) BETA - beta
  #
  # 5) SE - std. error
  #
  # 6) EAF – effect allele frequency
  #
  # 7) N - sample size
  #
  # 8) CHROMOSOME - chromosome of marker
  #
  # 9) POSITION - position of marker

  df_per_cohort <- per_cohort_summary_statistics_with_maf %>%
    select(MARKERNAME = variant,
           EA = eff_allele,
           NEA = non_eff_allele,
           BETA = beta,
           SE = standard_error,
           EAF = eaf,
           N = sample_size,
           CHROMOSOME = chromosome,
           POSITION = bp,
           cohort
    ) %>%
    group_by(cohort) %>%
    group_split()

  for (df in df_per_cohort) {
    cohort <- df$cohort[1]
    fwrite(df, sprintf("summary_stats_%s.txt.gz", cohort), col.names=T, row.names=F, quote=F)
  }

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}