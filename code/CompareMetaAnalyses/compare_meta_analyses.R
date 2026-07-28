#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)
library(argparse)
library(arrow)
library(metafor)


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
  help = "Path to resources.txt (tab-delimited key/value table)"
)

parser$add_argument(
  "--per_cohort_extractions", "-p",
  help = "RDS file with per cohort extractions")


# Define functions
meta_fit <- function(yi, vi) {
  res <- rma(yi, vi, method = "DL")
  tibble(
    beta = res$b[2],
    se = res$se,
    tau2 = res$tau2,
    I2 = res$I2
  )
}


forest_plot <- function(per_cohort_data, meta_analysed_data, hgnc_map, variant_reference, master_table, variant_gene_combinations, plot_legend = F) {

  per_cohort_results <- per_cohort_data %>% semi_join(
    variant_gene_combinations, by = c("phenotype" = "genes", "variant_index" = "variant_index")) %>% collect()
  meta_results <- meta_analysed_data %>% semi_join(
    variant_gene_combinations, by = c("phenotype" = "genes", "variant_index" = "variant_index")) %>% collect() %>%
      mutate(cohort = "meta")

  print("Extracting results done!")
  print(hgnc_map[variant_gene_combinations$genes])

  meta_and_per_cohort_bound <- bind_rows(per_cohort_results, meta_results) %>%
    inner_join(variant_reference) %>%
    mutate(symbol = hgnc_map[phenotype],
           symbol = factor(symbol, levels=unique(symbol), ordered=T)) %>%
    left_join(master_table) %>%
    mutate(
      sample_size = ifelse(!is.na(dataset_n), dataset_n, sample_size)) %>%
    arrange(desc(sample_size)) %>%
    mutate(
      cohort_new_name = replace_na(cohort_new_name, "Meta-analysis"),
      cohort_factor = factor(cohort, levels=rev(unique(cohort)), labels=rev(unique(cohort_new_name)), ordered=T),
      cohort_factor = fct_relevel(cohort_factor, "Meta-analysis", after = 0),
      rnaseq = case_when(cohort == "meta" ~ "Meta-analysis",
                         dataset_platform == "RNAseq" ~ "RNA-seq",
                         TRUE ~ "Microarray"),
      ancestry = ifelse(cohort == "meta", "Meta-analysis", ancestry)
    )

  meta_and_per_cohort_bound %>% distinct(cohort_factor)

  p <- ggplot(meta_and_per_cohort_bound, aes(x = beta, y = cohort_factor)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = (0.5 / (ggplot2::.pt * 72.27/96))) +
    geom_errorbarh(aes(xmin = beta-standard_error, xmax = beta+standard_error), height = 0.2, color = "gray40") +
    geom_point(size = 1.64, aes(color = ancestry, shape = rnaseq)) +
    scale_color_manual(values = pop_colors, na.value="black") +
    scale_shape_manual(values = c("Microarray" = 15, "Meta-analysis" = 23, "RNA-seq" = 17)) +
    facet_grid(~ variant + symbol, scales = "free_x") +
    labs(
      x = "Beta (+/- SE)",
      y = NULL,
      title = sprintf("Forest plots")
    ) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.line.x = element_line(color = "grey30", linewidth = (0.5 / (ggplot2::.pt * 72.27/96))),
      axis.ticks.x = element_line(color = "grey30", linewidth = (0.5 / (ggplot2::.pt * 72.27/96))),
      strip.background = element_rect(fill = "#f0f0f0", color = NA),
      strip.text = element_text(face = "bold"),
      aspect.ratio = 8,
      legend.position = if_else(plot_legend, "right", "none"),
      panel.border = element_rect(color = "grey10", linewidth = (0.5 / (ggplot2::.pt * 72.27/96)), fill=NA)
    )

  return(p)
}


generate_forest_plots <- function(
  per_cohort_data, meta_analysed_data, hgnc_map, variant_reference, master_table, example_variant_gene_combinations, grouping_var=c("variant_index")) {
  return(lapply(example_variant_gene_combinations %>% group_by(across(all_of(grouping_var))) %>% group_split(), function(group) {
    print(group)
    forest_plot(per_cohort_data, meta_analysed_data, hgnc_map,
                variant_reference = variant_reference, master_table = master_table,
                variant_gene_combinations = group, plot_legend = T)
  }))
}



# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  args <- parser$parse_args(argv)

  # Process input
  resources <- args$resources
  per_cohort_extractions <- args$per_cohort_extractions
  variables <- deframe(read_delim(resources, delim="\t"))

  mrmega_path <- "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/alternative_meta_analyses/mrmega_out_existing_genetic_axes.result"

  hgnc_genes <- fread("/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/network_analyses/hgnc_complete_set.txt")
  hgnc_map <- deframe(hgnc_genes %>% select(ensembl_gene_id, symbol))

  variant_reference <- read_parquet(variables["variant_reference"])

  cohort_info <- fread("/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/cohort_data/cohort_info.txt")
  master_table <- fread(variables["master_table"]) %>% select(cohort, cohort_new_name, dataset_n, dataset_platform, dataset_ancestry) %>%
    inner_join(cohort_info, by = c("cohort_new_name" = "cohort"))

  per_cohort_summary_statistics_raw <- readRDS(per_cohort_extractions)
  per_cohort_summary_statistics <- per_cohort_summary_statistics_raw %>%
    inner_join(variant_reference, by = "variant_index")

  # Load the main file of independent variants
  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type)) %>%
    select(variant_index, phenotype, beta, standard_error, type, typeCisTrans, sample_size, SusieRss_CS) %>%
    mutate(z_score = beta / standard_error,
           p_value = 2*pnorm(abs(z_score),lower.tail = F))

  signals_only_significant_in_one_cohort <- per_cohort_summary_statistics %>%
    mutate(p = 2*pnorm(abs(beta/standard_error), lower.tail = F)) %>%
    filter(!is.na(beta)) %>%
    group_by(variant_index, variant, phenotype) %>%
    summarise(n_cohorts = n(),
              n_cohorts_no_evidence = sum(p > 0.05),
              n_cohorts_evidence = n_cohorts - n_cohorts_no_evidence,
              cohorts_significant = paste0(cohort[p <= 0.05], collapse=";"),
              significance_in_one_cohort = min(p)) %>%
    filter(n_cohorts_evidence < 2) %>%
    inner_join(independent_variant_set_annot_filtered, by = join_by(variant_index, phenotype)) %>%
    filter(!is.na(SusieRss_CS))

  signals_only_significant_in_one_cohort %>% group_by(cohorts_significant) %>% tally()

  mrmega <- fread(mrmega_path) %>%
    separate_wider_delim(MarkerName, "__", names = c("variant", "phenotype")) %>%
    inner_join(variant_reference, by = c("variant"))

  mrmega_joined <- mrmega %>%
    inner_join(independent_variant_set_annot_filtered, by =c("variant_index", "phenotype")) %>%
    select(-`P-value_ancestry_het`, -`P-value_residual_het`, -`P-value_association`) %>%
    mutate(mlog10p_fixed = -log10(pmax(p_value, .Machine$double.xmin)),
           z_score_abs = abs(z_score),
           z_score_mrmega = sqrt(chisq_association),
           z_diff = z_score_mrmega - z_score_abs,
           p_value_pt = 2*pt(z_score_abs, sample_size - 1, lower.tail = FALSE, log.p = F),
           p_value_association = pchisq(chisq_association, ndf_association, lower.tail = FALSE, log.p = F),
           mlog10p_mrmega = -log10(pmax(p_value_association, .Machine$double.xmin)),
           pvalue_residual_het = pchisq(chisq_residual_het, ndf_residual_het, lower.tail = FALSE, log.p = F),
           pvalue_ancestry_het = pchisq(chisq_ancestry_het, ndf_ancestry_het, lower.tail = FALSE, log.p = F)
    )

  mrmega_joined %>%
    mutate(largest_heterogeneity_factor = if_else(chisq_residual_het > chisq_ancestry_het, "Residual", "Ancestry")) %>%
    group_by(largest_heterogeneity_factor, typeCisTrans) %>%
    tally() %>%
    group_by(typeCisTrans) %>%
    mutate(p = n / sum(n))

  ggplot(mrmega_joined, aes(chisq_residual_het, chisq_ancestry_het)) +
    geom_abline() +
    geom_smooth(method="lm") +
    geom_point(shape=16, alpha=0.1, size=(1 / (ggplot2::.pt * 72.27 / 96))) +
    coord_fixed() +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(vars(typeCisTrans))

  ggsave("driving_ancestry.pdf")

  fwrite(mrmega, "mrmega_20260623.txt.gz", sep="\t", row.names=F, col.names=T)

  correlation <- cor(mrmega_joined$z_score_abs, mrmega_joined$z_score_mrmega)

  mrmega_joined %>%
    mutate(z_sig = (z_diff / sqrt(2 - 2*correlation))) %>%
    select(variant,phenotype,z_score_abs, z_score_mrmega, z_diff, z_sig, p_values_fixed, p_values_mrmega) %>%
    summarise(n_more_significant = sum(z_diff > 0),
              n_less_significant = sum(z_diff < 0),
              )

  ggplot(mrmega_joined, aes(z_score^2, chisq_association)) +
    xlab("Chi2 statistic (Fixed-effects meta-analysis)") +
    ylab("Chi2 statistic (MR-MEGA meta-analysis)") +
    geom_point(alpha=0.2, shape=16, size=1)

  ggsave("mrmega_vs_fixed_chi2.pdf", height=4, width=4)

  ggplot(mrmega_joined, aes(mlog10p_fixed, mlog10p_mrmega2)) +
    xlab("-log10(p) (Fixed-effects meta-analysis)") +
    ylab("-log10(p) (MR-MEGA meta-analysis)") +
    geom_point(alpha=0.2, shape=16, size=1)

  ggsave("mrmega_joined.pdf", height=4, width=4)

  mrmega_joined %>% group_by(pvalue_ancestry_het > 0.0001, typeCisTrans) %>%
    tally() %>%
    group_by(typeCisTrans) %>%
    mutate(p = n/sum(n)*100)

  ggplot(mrmega_joined, aes(typeCisTrans, -log10(pmax(pvalue_ancestry_het, 1e-16)))) +
    geom_violin()

  ggsave("mrmega_ancestry_het.pdf", height=4, width=4)

  mrmega_joined %>% group_by(pvalue_residual_het > 0.0001, typeCisTrans) %>%
    tally() %>%
    group_by(typeCisTrans) %>%
    mutate(p = n/sum(n)*100)

  ggplot(mrmega_joined, aes(typeCisTrans, -log10(pmax(pvalue_residual_het, 1e-16)))) +
    geom_violin()

  ggsave("mrmega_residual_het.pdf", height=4, width=4)

  mrmega_combinations <- mrmega_joined %>%
    filter(pvalue_residual_het == 0, typeCisTrans == "trans") %>%
    select(variant_index, genes = phenotype) %>% slice_head(n=5) %>% mutate(page=1)

  residual_het_trans_forest_plots <- generate_forest_plots(
    per_cohort_summary_statistics_raw, independent_variant_set_annot_filtered,
    hgnc_map, variant_reference, master_table, mrmega_combinations, grouping_var = c("page"))

  pdf("residual_het_trans_forest_plots.pdf", onefile = TRUE, width=180/10/2.54, height = 180/10/2.54)
  for (i in seq(length(residual_het_trans_forest_plots))) {
    print(residual_het_trans_forest_plots[[i]])
  }
  dev.off()

  mrmega_combinations_ancestry <- mrmega_joined %>%
    filter(typeCisTrans == "trans") %>%
    arrange(desc(chisq_ancestry_het)) %>%
    slice_head(n=20) %>%
    mutate(page = (row_number() - 1) %/% 4) %>%
    select(variant_index, genes = phenotype, pvalue_ancestry_het, page)

  ancestry_het_trans_forest_plots <- generate_forest_plots(
    per_cohort_summary_statistics_raw, independent_variant_set_annot_filtered,
    hgnc_map, variant_reference, master_table, mrmega_combinations_ancestry, grouping_var = "page")

  pdf("ancestry_het_trans_forest_plots.pdf", onefile = TRUE, width=180/10/2.54, height = 180/10/2.54)
  for (i in seq(length(ancestry_het_trans_forest_plots))) {
    print(ancestry_het_trans_forest_plots[[i]])
  }
  dev.off()

  # Work
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}