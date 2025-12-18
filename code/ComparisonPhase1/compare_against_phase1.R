#!/usr/bin/env Rscript


# Load libraries
library(data.table)
library(tidyverse)
library(arrow)
library(broom)
library(dbplyr)
library(ggrastr)

# Declare constants
old <- theme_set(theme_classic(base_size = 8))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27 / 96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  text = element_text(family = "Helvetica"),
  axis.title = element_text(colour = "black"),
  title = element_text(colour = "#595A5C", face = "plain")
)


# Declare function definitions
#' Function for Rb analysis
#'
#' @param b1 Beta from first dataset.
#' @param se1 Standard error of beta from first dataset.
#' @param b2 Beta from second dataset.
#' @param se2 Standard error of beta from second dataset.
#' @param theta Variable representing sample overlap between two datasets. Should be set 0 if no sample overlap.
#'
#' @return Data frame with Rb, SE(Rb) and corresponding P-value.
#' @export
#'
#' @note This function is slightly adapted from the script shared by Ting Qi.
#'
#' @examples
calcu_cor_true <- function(b1, se1, b2, se2, theta) {
  idx <- which(is.infinite(b1) |
                 is.infinite(b2) |
                 is.infinite(se1) |
                 is.infinite(se2))
  if (length(idx) > 0) {
    b1 <- b1[-idx]
    se1 <- se1[-idx]
    b2 <- b2[-idx]
    se2 <- se2[-idx]
    theta <- theta[-idx]
  }

  var_b1 <- var(b1, na.rm = T) - mean(se1^2, na.rm = T)
  var_b2 <- var(b2, na.rm = T) - mean(se2^2, na.rm = T)
  if (var_b1 < 0) {
    var_b1 <- var(b1, na.rm = T)
  }
  if (var_b2 < 0) {
    var_b2 <- var(b2, na.rm = T)
  }
  cov_b1_b2 <- cov(b1, b2, use = "complete.obs") - mean(theta, na.rm = T) * sqrt(mean(se1^2, na.rm = T) * mean(se2^2, na.rm = T))
  r <- cov_b1_b2 / sqrt(var_b1 * var_b2)

  r_jack <- c()
  n <- length(b1)
  for (k in 1:n) {
    b1_jack <- b1[-k]
    se1_jack <- se1[-k]
    var_b1_jack <- var(b1_jack, na.rm = T) - mean(se1_jack^2, na.rm = T)
    b2_jack <- b2[-k]
    se2_jack <- se2[-k]
    var_b2_jack <- var(b2_jack, na.rm = T) - mean(se2_jack^2, na.rm = T)
    if (var_b1_jack < 0) {
      var_b1_jack <- var(b1_jack, na.rm = T)
    }
    if (var_b2_jack < 0) {
      var_b2_jack <- var(b2_jack, na.rm = T)
    }
    theta_jack <- theta[-k]
    cov_e1_jack_e2_jack <- mean(theta_jack, na.rm = T) * sqrt(mean(se1_jack^2, na.rm = T) * mean(se2_jack^2, na.rm = T))
    cov_b1_b2_jack <- cov(b1_jack, b2_jack, use = "complete.obs") - cov_e1_jack_e2_jack
    r_tmp <- cov_b1_b2_jack / sqrt(var_b1_jack * var_b2_jack)
    r_jack <- c(r_jack, r_tmp)
  }
  r_mean <- mean(r_jack, na.rm = T)
  idx <- which(is.na(r_jack))
  if (length(idx) > 0) {
    se_r <- sqrt((n - 1) / n * sum((r_jack[-idx] - r_mean)^2))
  } else {
    se_r <- sqrt((n - 1) / n * sum((r_jack - r_mean)^2))
  }

  p <- pchisq((r / se_r)**2, df = 1, lower.tail = FALSE)

  res <- data.frame("r" = r, "se_r" = se_r, "p" = p)
  return(res)
}

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))

  # Load variant annotations, has various annotations besides just positions, names and alleles
  variant_reference <- read_parquet(variables["variant_reference"])

  # Process input
  meta_analysed <- fread(variables["independent_variant_set"]) %>%
    select(phenotype, variant_index, beta, standard_error, type, p_value)

  phase2_ds <- open_dataset("") # Here goes the folder with the eQTLGen meta-analysis output

  phase1_ds <- open_dataset("harmonized/eqtlgen_parquet_dataset") %>%
    select(Gene, variant_index, beta, standard_error = se, type = gwas_id, p)

  phase1_sign_ds <- open_dataset("harmonized_fdr_pass/eqtlgen_parquet_dataset") %>%
    select(Gene, variant_index, beta, standard_error = se, type = gwas_id, p)

  max_p_cis <- phase1_sign_ds %>%
    filter(type == "eqtlgen_phase1_cis") %>%
    collect() %>%
    pull(p) %>%
    max()
  max_p_trans <- phase1_sign_ds %>%
    filter(type == "eqtlgen_phase1_trans") %>%
    collect() %>%
    pull(p) %>%
    max()

  lead_cis_eQTL_effects <- phase1_ds %>%
    filter(type == "eqtlgen_phase1_cis") %>%
    collect() %>%
    group_by(Gene) %>%
    slice_min(p, with_ties=F) %>%
    mutate(FDR_significant = p <= max_p_cis) %>%
    inner_join(variant_reference)

  fwrite(lead_cis_eQTL_effects, "lead_cis-eQTL_per_gene_eQTLGen_phase1_harmonized_to_phase2.txt.gz",
         sep = "\t", col.names = T, row.names = F, quote = F)

  minimum_p_per_gene_cis <- phase1_ds %>%
    filter(type == "eqtlgen_phase1_cis") %>%
    group_by(variant_index) %>%
    summarise(minimum_p_value_over_genes = min(p)) %>%
    collect() %>%
    mutate(any_fdr_significant = minimum_p_value_over_genes <= max_p_cis) %>%
    inner_join(variant_reference)

  fwrite(minimum_p_per_gene_cis, "minimum_p_value_per_gene_eQTLGen_phase1_harmonized_to_phase2.txt.gz",
         sep = "\t", col.names = T, row.names = F, quote = F)

  phase1_joined <- phase1_ds %>%
    right_join(
      meta_analysed,
      by = c("Gene" = "phenotype", "variant_index" = "variant_index"),
      suffix = c("_phase1", "_phase2")) %>%
    collect() %>%
    mutate(
      z_phase1 = beta_phase1 / standard_error_phase1,
      z_phase2 = beta_phase2 / standard_error_phase2,
    )

  phase1_joined_filtered <- phase1_joined %>% filter(!is.na(type_phase1))

  phase1_joined_filtered %>%
    group_by(type_phase1) %>%
    summarise(n = n(), alpha = 0.05 / n, significant_for_number_of_tests = sum(p < alpha), nominal = sum(p < 0.05), p_bonferroni = significant_for_number_of_tests / n * 100, p_nominal = nominal / n * 100)

  phase1_joined_filtered %>%
    filter(p_value < 5e-8) %>%
    group_by(type_phase1) %>%
    summarise(n = n(), alpha = 0.05 / n, significant_for_number_of_tests = sum(p < alpha), nominal = sum(p < 0.05), p_bonferroni = significant_for_number_of_tests / n * 100, p_nominal = nominal / n * 100)

  cor_stats <- phase1_joined_filtered %>%
    group_by(type_phase1) %>%
    summarise(tidy(cor.test(beta_phase1, beta_phase2, use = "pairwise.complete.obs")))


  #rb_mat = calcu_cor_true(beta_phase1, standard_error_phase1, beta_phase2, standard_error_phase2, theta = 1),

  cor_stats_rb <- phase1_joined_filtered %>%
    group_by(type_phase1) %>%
    summarise(rb_mat = calcu_cor_true(beta_phase1, standard_error_phase1, beta_phase2, standard_error_phase2, theta = 1),
              z_r = tidy(cor.test(z_phase1, z_phase2))) %>%
    unnest(cols = c(rb_mat, z_r), names_sep="_") %>%
    mutate(rb_lab = sprintf("Rb of betas: %.3f (SE = %.3f)", rb_mat_r, rb_mat_se_r),
           z_r_lab = sprintf("R of Z-scores: %.3f (SE = %.3f)", z_r_estimate, z_r_estimate - z_r_conf.low),
           annot_lab = paste0(rb_lab, "\n", z_r_lab))

  cor_stats_rb_trans_type <- phase1_joined_filtered %>%
    filter(type_phase1 == "eqtlgen_phase1_trans") %>%
    group_by(ctc_eqtl) %>%
    summarise(
      nest(calcu_cor_true(beta_phase1, standard_error_phase1, beta_phase2, standard_error_phase2, theta = 1)), n = n()) %>%
    unnest_wider(data) %>%
    mutate(rb_lab = sprintf("%.2f (%.2f, %.2e)", r, se_r, p))

  colnames(cor_stats_rb) <- c("type_phase1", "estimate", "se", "p")

  phase1_joined_filtered_cis <- phase1_joined_filtered %>% filter(type_phase1 == "eqtlgen_phase1_cis")

  calcu_cor_true(phase1_joined_filtered_cis$beta_phase1, phase1_joined_filtered_cis$standard_error_phase1,
                 phase1_joined_filtered_cis$beta_phase2, phase1_joined_filtered_cis$standard_error_phase2, theta = 10000)

  cor_plots <- ggplot(phase1_joined %>% filter(!is.na(type_phase1)), aes(beta_phase1, beta_phase2)) +
    rasterize(geom_point(alpha = 0.1, shape=16), dpi=300) +
    geom_abline() +
    geom_smooth(method = "lm") +
    labs(title = "Correlation of betas") +
    xlab("Phase 1 Betas") +
    ylab("Phase 2 Betas") +
    geom_text(data = cor_stats, aes(x = -Inf, y = Inf, label = estimate), vjust = "inward", hjust = "inward") +
    facet_wrap(vars(type_phase1), scales = "free") +
    theme(aspect.ratio = 1)

  ggsave("correlations_to_phase1_betas.png", cor_plots)

  cor_plots <- ggplot(
    phase1_joined_filtered %>% filter(type_phase1 == "eqtlgen_phase1_trans"),
    aes(beta_phase1, beta_phase2, color = ctc_eqtl)) +
    rasterize(geom_point(alpha = 0.1, shape=16), dpi=300) +
    geom_abline() +
    geom_smooth(method = "lm") +
    labs(title = "Correlation of betas") +
    xlab("Phase 1 Betas") +
    ylab("Phase 2 Betas") +
    scale_color_manual(values = c("#D55E00", "black")) +
    geom_text(data = cor_stats_rb_trans_type %>%
      ungroup() %>%
      summarise(rb_lab = paste0(ctc_eqtl, " | ", rb_lab, collapse = "\n")), aes(x = -Inf, y = Inf, label = rb_lab), vjust = "inward", hjust = "inward", inherit.aes = F) +
    facet_wrap(vars(ctc_eqtl), scales = "free") +
    theme(aspect.ratio = 1)

  ggsave("correlations_to_phase1_betas_split_trans.png", cor_plots)


  cor_stats_z <- phase1_joined %>%
    filter(!is.na(type_phase1)) %>%
    group_by(type_phase1) %>%
    summarise(tidy(cor.test(z_phase1, z_phase2, use = "pairwise.complete.obs")))

  cor_plots_z <- ggplot(phase1_joined %>% filter(!is.na(type_phase1)), aes(z_phase1, z_phase2)) +
    rasterize(geom_point(alpha = 0.1, shape=16), dpi=300) +
    geom_abline() +
    geom_smooth(method = "lm") +
    labs(title = "Correlation of Z-scores") +
    xlab("Phase 1 Z-scores") +
    ylab("Phase 2 Z-scores") +
    geom_text(data = cor_stats_rb, aes(x = -Inf, y = Inf, label = annot_lab), vjust = "inward", hjust = "inward") +
    facet_wrap(vars(type_phase1), scales = "free") +
    theme(aspect.ratio = 1)

  ggsave("scatterplot_lookup_phase2_in_phase1_zscores.pdf", cor_plots_z)

  # Now, focus on the opposite comparison, in which phase 1 is discovery
  phase1_trans_fdr_pass_df <- phase1_sign_ds %>%
    filter(type == "eqtlgen_phase1_trans") %>%
    collect() %>%
    mutate(variant_index = as.integer64(variant_index))

  phase1_trans_fdr_pass_joined <- phase2_ds %>%
    right_join(phase1_trans_fdr_pass_df, by = c("phenotype" = "Gene", "variant_index" = "variant_index")) %>%
    collect()


  phase1_trans_fdr_pass_df <- phase1_sign_ds %>%
    filter(type == "eqtlgen_phase1_trans") %>%
    collect()
  phase1_trans_fdr_pass_df_lookup <- phase1_trans_fdr_pass_df %>%
    group_by(Gene) %>%
    summarise(variant_indices = list(variant_index)) %>%
    mutate(row_num = row_number())

  # Per failed gene-window combination, get the lead variant.
  queried_phase1_eqtls <- pmap(
    phase1_trans_fdr_pass_df_lookup,
    function(row_num, Gene, variant_indices) {
      print(sprintf("%s: %s, %s", row_num, Gene, length(variant_indices)))
      phase2_ds %>%
        filter(
          phenotype == !!Gene,
          variant_index %in% !!variant_indices
        ) %>%
        collect()
    }
  )

  queried_phase1_trans_eqtls_bound <- bind_rows(queried_phase1_eqtls)
  saveRDS(queried_phase1_trans_eqtls_bound, "phase1_trans_eqlts_queried_from_phase2.rds")

  phase1_cis_fdr_pass_df <- phase1_sign_ds %>%
    filter(type == "eqtlgen_phase1_cis") %>%
    collect() %>%
    group_by(Gene) %>%
    slice_max(abs(beta/standard_error), with_ties = F) %>%
    mutate(row_num = row_number())

  # Per failed gene-window combination, get the lead variant.
  queried_phase1_eqtls_cis <- pmap(
    phase1_cis_fdr_pass_df %>% select(row_num, Gene, variant_index),
    function(row_num, Gene, variant_index) {
      print(sprintf("%s: %s, %s", row_num, Gene, length(variant_index)))
      phase2_ds %>%
        filter(
          phenotype == !!Gene,
          variant_index %in% !!variant_index
        ) %>%
        collect()
    }
  )

  queried_phase1_cis_eqtls_bound <- bind_rows(queried_phase1_eqtls_cis)
  saveRDS(queried_phase1_cis_eqtls_bound, "top_phase1_cis_eqlts_queried_from_phase2_20251119.rds")


  queried_phase1_trans_eqtls_bound <- readRDS("phase1_trans_eqlts_queried_from_phase2.rds")
  queried_phase1_cis_eqtls_bound <- readRDS("top_phase1_cis_eqlts_queried_from_phase2_20251119.rds")

  phase1_cis_fdr_pass_df_with_phase2 <- phase1_cis_fdr_pass_df %>%
    left_join(
      queried_phase1_cis_eqtls_bound,
      by = c("Gene" = "phenotype", "variant_index" = "variant_index"), suffix=c("_phase1", "_phase2")) %>%
    mutate(z_score_phase2 = beta_phase2 / standard_error_phase2, p_value_phase2 = 2 * pnorm(abs(z_score_phase2), lower.tail = F),
           z_score_phase1 = beta_phase1 / standard_error_phase1)

  phase1_cis_fdr_pass_df_with_phase2 %>%
    summarise(
      is_available = sum(!is.na(p_value)),
      is_gw_significant = sum(p_value < 5e-8, na.rm = T),
      is_p1_significant = sum(p_value < max_p_cis, na.rm = T),
      is_nom_significant = sum(p_value < 0.05, na.rm = T))

  phase1_trans_fdr_pass_df_with_phase2 <- phase1_trans_fdr_pass_df %>%
    left_join(
      queried_phase1_trans_eqtls_bound,
      by = c("Gene" = "phenotype", "variant_index" = "variant_index"), suffix=c("_phase1", "_phase2")) %>%
    mutate(z_score_phase2 = beta_phase2 / standard_error_phase2, p_value_phase2 = 2 * pnorm(abs(z_score_phase2), lower.tail = F),
           z_score_phase1 = beta_phase1 / standard_error_phase1)

  phase1_trans_fdr_pass_df_with_phase2 %>%
    summarise(
      is_available = sum(!is.na(p_value)),
      is_gw_significant = sum(p_value < 5e-8, na.rm = T),
      is_p1_significant = sum(p_value < max_p_trans, na.rm = T),
      is_nom_significant = sum(p_value < 0.05, na.rm = T))

  phase1_lookup_in_phase2_joined <- bind_rows(list("cis" = phase1_cis_fdr_pass_df_with_phase2, "trans"=phase1_trans_fdr_pass_df_with_phase2), .id="type_phase1")

  cor_stats_rb_reverse <- phase1_lookup_in_phase2_joined %>%
    group_by(type_phase1) %>%
    summarise(rb_mat = calcu_cor_true(beta_phase1, standard_error_phase1, beta_phase2, standard_error_phase2, theta = 1),
              z_r = tidy(cor.test(z_score_phase1, z_score_phase2)),
              n = sum(!is.na(beta_phase1) & !is.na(beta_phase2))) %>%
    unnest(cols = c(rb_mat, z_r), names_sep="_") %>%
    mutate(rb_lab = sprintf("Rb of betas: %.3f (SE = %.3f)", rb_mat_r, rb_mat_se_r),
           z_r_lab = sprintf("R of Z-scores: %.3f (SE = %.3f)", z_r_estimate, z_r_estimate - z_r_conf.low),
           annot_lab = paste0(sprintf("n = %s", n), "\n", rb_lab, "\n", z_r_lab))

  cor_plots_b_lookup_phase1_in_phase2 <- ggplot(phase1_lookup_in_phase2_joined %>% filter(!is.na(type_phase1)), aes(beta_phase1, beta_phase2)) +
    rasterize(geom_point(alpha = 0.1, shape=16), dpi=300) +
    geom_abline() +
    geom_smooth(method = "lm") +
    labs(title = "Correlation of betas") +
    xlab("Phase 1 Betas") +
    ylab("Phase 2 Betas") +
    geom_text(data = cor_stats_rb_reverse, aes(x = -Inf, y = Inf, label = rb_lab), vjust = "inward", hjust = "inward") +
    facet_wrap(vars(type_phase1), scales = "free") +
    theme(aspect.ratio = 1)

  ggsave("scatterplot_phase1_lookup_in_phase2_betas.pdf", cor_plots_b_lookup_phase1_in_phase2)


  cor_plots_z_lookup_phase1_in_phase2 <- ggplot(phase1_lookup_in_phase2_joined %>% filter(!is.na(type_phase1)), aes(z_score_phase1, z_score_phase2)) +
    rasterize(geom_point(alpha = 0.1, shape=16), dpi=300) +
    geom_abline() +
    geom_smooth(method = "lm") +
    labs(title = "Correlation of Z-scores") +
    xlab("Phase 1 Z-scores") +
    ylab("Phase 2 Z-scores") +
    geom_text(data = cor_stats_rb_reverse, aes(x = -Inf, y = Inf, label = annot_lab), vjust = "inward", hjust = "inward") +
    facet_wrap(vars(type_phase1), scales = "free") +
    theme(aspect.ratio = 1)

  ggsave("scatterplot_phase1_lookup_in_phase2_zscores_updated.pdf", cor_plots_z_lookup_phase1_in_phase2)


  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}