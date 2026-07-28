#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)
library(arrow)
library(broom)
library(tidyr)


# Declare constants

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

  # Work
  master_table <- fread(variables["master_table"])

  master_table_eur <- master_table %>% filter(str_detect(dataset_ancestry, "EUR|Askenazi"))
  master_table_non_eur <- master_table %>% filter(!str_detect(dataset_ancestry, "EUR|Askenazi"))



  values <- master_table_eur$dataset_n
  names(values) <- master_table_eur$cohort_new_name
  target <- 31684
  target2 <- 35988

# function to evaluate combinations of size k
  find_combos <- function(k, target) {
    cmb <- combn(seq_along(values), k)
    sums <- colSums(matrix(values[cmb], nrow = k))
    diff <- target - sums

    data.frame(
      datasets = apply(cmb, 2, function(i) {paste(names(values[i]), collapse=",")}),
      rows = apply(cmb, 2, paste, collapse = ","),
      sum = sums,
      diff = diff
    )
  }

  res <- bind_rows(lapply(37, find_combos, target=target))

  res_filtered <- res[res$diff == 0,]

  res_long <- res_filtered %>% mutate(combination_nr = row_number()) %>%
    separate_longer_delim(c(datasets, rows), ",") %>%
    inner_join(master_table_eur, by = c("datasets" = "cohort_new_name"))

  res_summarised <- res_long %>%
    group_by(combination_nr) %>%
    summarise(n_cohorts=n(), n_samples=sum(dataset_n), n_p1 = sum(phase1), p_p1 = n_p1/n_cohorts)

  datasets_eur_subsampled_to_p1 <- res_long %>% filter(combination_nr == 47) %>% pull(datasets)
  datasets_non_eur <- master_table_non_eur %>% pull(cohort_new_name)
  datasets_eur <- master_table_eur %>% pull(cohort_new_name)
  datasets_eur_subsampled_to_p1_and_non_eur <- c(datasets_eur_subsampled_to_p1, datasets_non_eur)
  dput(datasets_eur_subsampled_to_p1)
  dput(datasets_eur)

  res <- bind_rows(lapply(39, find_combos, target=target2))
  print(res[which.min(abs(res$diff)),])

  res_filtered <- res[which.min(abs(res$diff)),]

  res_long <- res_filtered %>% mutate(combination_nr = row_number()) %>%
    separate_longer_delim(c(datasets, rows), ",") %>%
    inner_join(master_table_eur, by = c("datasets" = "cohort_new_name"))

  datasets_eur_subsampled_to_multi <- res_long %>% pull(datasets)
  dput(datasets_eur_subsampled_to_multi)

  datasets_eur_subsampled_to_multi <- c("DIRECT_NonFIN", "DIRECT_FIN", "Cruchagalab_EUR", "KORA_batch2",
                                        "EstBB_IlluminaHT12v3", "EstBB_RNAseq", "KORA_batch1", "GTEx_v9_EUR",
                                        "YoungFinns", "CHDWB", "CAD", "OphoffBP_GSA_lowcov", "OphoffBP_OMEX_highcov",
                                        "OphoffBP_OMEX_lowcov", "OphoffBP_PSYCH_lowcov", "GAIT2", "PRECISESADS",
                                        "GAinS_RNAseq", "SHIP_TREND", "LIFE_HEART", "DGN", "STRIP", "CanPath_freeze1",
                                        "CanPath_freeze2", "AMP-PD_AJ", "AMP-PD_EUR", "NESDA_NTR", "Rotterdam_IlluminaHT12v4",
                                        "Fehrmann-H8v2", "Fehrmann-HT12v3", "CODAM", "LL", "LLS_660Q",
                                        "LLS_OmniExpr", "NTR_Affy6", "NTR_GoNL", "PAN", "Rotterdam_RNAseq",
                                        "BSGS", "INTERVAL", "InCHIANTI")

  datasets_eur_subsampled_to_p1 <- c("DIRECT_NonFIN", "Cruchagalab_EUR", "KORA_batch2", "EstBB_IlluminaHT12v3",
                                     "EstBB_RNAseq", "KORA_batch1", "GTEx_v9_EUR", "YoungFinns", "CHDWB",
                                     "CAD", "OphoffBP_GSA_lowcov", "OphoffBP_OMEX_lowcov", "OphoffBP_PSYCH_lowcov",
                                     "PRECISESADS", "SHIP_TREND", "LIFE_ADULT", "LIFE_HEART", "DGN",
                                     "STRIP", "CanPath_freeze1", "CanPath_freeze2", "AMP-PD_AJ", "AMP-PD_EUR",
                                     "NESDA_NTR", "Rotterdam_IlluminaHT12v4", "Fehrmann-H8v2", "Fehrmann-HT12v3",
                                     "CODAM", "LL", "LLS_660Q", "LLS_OmniExpr", "NTR_Affy6", "NTR_GoNL",
                                     "PAN", "Rotterdam_RNAseq", "BSGS", "InCHIANTI")

  datasets_eur <- c("DIRECT_NonFIN", "DIRECT_FIN", "Cruchagalab_EUR", "KORA_batch2",
                    "EstBB_IlluminaHT12v3", "EstBB_RNAseq", "KORA_batch1", "GTEx_v9_EUR",
                    "YoungFinns", "CHDWB", "CAD", "OphoffBP_GSA_lowcov", "OphoffBP_OMEX_highcov",
                    "OphoffBP_OMEX_lowcov", "OphoffBP_PSYCH_lowcov", "GAIT2", "PRECISESADS",
                    "GAinS_RNAseq", "SHIP_TREND", "LIFE_ADULT", "LIFE_HEART", "DGN",
                    "STRIP", "CanPath_freeze1", "CanPath_freeze2", "AMP-PD_AJ", "AMP-PD_EUR",
                    "NESDA_NTR", "Rotterdam_IlluminaHT12v4", "Fehrmann-H8v2", "Fehrmann-HT12v3",
                    "CODAM", "LL", "LLS_660Q", "LLS_OmniExpr", "NTR_Affy6", "NTR_GoNL",
                    "PAN", "Rotterdam_RNAseq", "BSGS", "INTERVAL", "InCHIANTI")

  datasets_non_eur <- master_table_non_eur %>% pull(cohort_new_name)
  datasets_eur <- master_table_eur %>% pull(cohort_new_name)
  datasets_eur_subsampled_to_p1_and_non_eur <- c(datasets_eur_subsampled_to_p1, datasets_non_eur)

  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type)) %>%
    select(variant_index, phenotype, beta, standard_error, type, typeCisTrans, locus_str) %>%
    mutate(z_score = beta / standard_error) %>%
    group_by(phenotype, locus_str) %>%
    slice_max(abs(z_score), with_ties=F)

  per_cohort_extractions <- variables["per_cohort_extractions"]
  per_cohort_summary_statistics <- readRDS(per_cohort_extractions) %>%
    semi_join(independent_variant_set_annot_filtered, by = c("variant_index", "phenotype"))
  per_cohort_summary_statistics <- per_cohort_summary_statistics %>%
    mutate(z_score = beta / standard_error) %>%
    left_join(master_table %>% select(cohort, cohort_new_name), by = c("cohort" = "cohort"))

  meta_analysis_eur_subsampled <- per_cohort_summary_statistics %>%
    filter(cohort_new_name %in% datasets_eur_subsampled_to_p1, !is.na(beta)) %>%
    group_by(phenotype, variant_index) %>%
    mutate(weights = 1 / standard_error^2) %>%
    summarise(
      eur_31k_beta = sum(beta * weights) / sum(weights),
      eur_31k_se = sqrt(1 / sum(weights)),
      eur_31k_z = eur_31k_beta / eur_31k_se
    )

  meta_analysis_eur_all <- per_cohort_summary_statistics %>%
    filter(cohort_new_name %in% datasets_eur, !is.na(beta)) %>%
    group_by(phenotype, variant_index) %>%
    mutate(weights = 1 / standard_error^2) %>%
    summarise(
      eur_all_beta = sum(beta * weights) / sum(weights),
      eur_all_se = sqrt(1 / sum(weights)),
      eur_all_z = eur_all_beta / eur_all_se
    )

  meta_analysis_subsampled_and_non_eur <- per_cohort_summary_statistics %>%
    filter(cohort_new_name %in% datasets_eur_subsampled_to_p1_and_non_eur, !is.na(beta)) %>%
    group_by(phenotype, variant_index) %>%
    mutate(weights = 1 / standard_error^2) %>%
    summarise(
      multi_beta = sum(beta * weights) / sum(weights),
      multi_se = sqrt(1 / sum(weights)),
      multi_z = multi_beta / multi_se
    )

  meta_analysis_35k_subsampled <- per_cohort_summary_statistics %>%
    filter(cohort_new_name %in% datasets_eur_subsampled_to_multi, !is.na(beta)) %>%
    group_by(phenotype, variant_index) %>%
    mutate(weights = 1 / standard_error^2) %>%
    summarise(
      eur_35k_beta = sum(beta * weights) / sum(weights),
      eur_35k_se = sqrt(1 / sum(weights)),
      eur_35k_z = eur_35k_beta / eur_35k_se
    )

  summary_statistics_merged <- independent_variant_set_annot_filtered %>%
    left_join(meta_analysis_eur_subsampled, by = c("variant_index", "phenotype")) %>%
    left_join(meta_analysis_eur_all, by = c("variant_index", "phenotype"))

  z_score_threshold <- abs(qnorm(5e-8/2))


  # We can quantify: the number of primary eQTLs found, the
  summary_statistics_diff <- summary_statistics_merged %>%
    ungroup() %>%
    summarise(
      n = n(),
      all_singificant = sum(abs(z_score) > z_score_threshold),
      eur_all_significant = sum(abs(eur_all_z) > z_score_threshold),
      eur_31k_significant = sum(abs(eur_31k_z) > z_score_threshold),
      eur_all_significant_p = all_singificant / eur_all_significant,
      eur_31k_significant_p = all_singificant / eur_31k_significant,
      average_increase_compared_to_eur_31k = mean(z_score / eur_31k_z),
      average_increase_compared_to_eur_all = mean(z_score / eur_all_z))

  summary_statistics_diff %>% pivot_longer(cols=everything())

  # We can quantify: the number of primary eQTLs found, stratified by cis-vs-trans
  summary_statistics_diff <- summary_statistics_merged %>%
    ungroup() %>%
    group_by(typeCisTrans) %>%
    summarise(
      n = n(),
      all_singificant = sum(abs(z_score) > z_score_threshold),
      eur_all_significant = sum(abs(eur_all_z) > z_score_threshold),
      eur_31k_significant = sum(abs(eur_31k_z) > z_score_threshold),
      eur_all_significant_p = all_singificant / eur_all_significant,
      eur_31k_significant_p = all_singificant / eur_31k_significant,
      average_increase_compared_to_eur_31k = mean(z_score / eur_31k_z),
      average_increase_compared_to_eur_all = mean(z_score / eur_all_z))

  summary_statistics_diff %>% pivot_longer(cols=-typeCisTrans)


  # From all eQTLGen phase 2 cohorts, we assessed all possible combinations of 37 cohorts that had a total sample size of 31684.
  # We arbitrarily picked one of the two cohort combinations that had the maximum number (26) of overlapping cohorts from eQTLGen phase 1.

  # For these, we extracted the per-cohort summary statistics for all top finemapped eQTLs per locus-gene combination and
  # meta-analysed the beta's and standard errors using this subset of cohorts.

  # Compared to a meta-analysis subsampled to 31K european-only samples, p2's primary analysis yields 1.55 times more significant primary eQTLs.
  # In the full meta-analysis, absolute Z-scores are approximately 1.22 times larger compared to the subsampled analysis of 31K european-only cohorts.
  # Compared to a meta-analysis subsampled to all european samples, p2's primary analysis yields 1.13 times more significant primary eQTLs.
  # In the full meta-analysis, absolute Z-scores are approximately 1.04 times larger compared to the subsampled analysis of only european cohorts.



  summary_statistics_merged_35k <- independent_variant_set_annot_filtered %>%
    left_join(meta_analysis_eur_subsampled, by = c("variant_index", "phenotype")) %>%
    left_join(meta_analysis_35k_subsampled, by = c("variant_index", "phenotype")) %>%
    left_join(meta_analysis_subsampled_and_non_eur, by = c("variant_index", "phenotype"))

  # We can quantify: the number of primary eQTLs found, the
  summary_statistics_diff <- summary_statistics_merged_35k %>%
    ungroup() %>%
    summarise(
      n = n(),
      all_singificant = sum(abs(z_score) > z_score_threshold),
      multi_significant = sum(abs(multi_z) > z_score_threshold),
      eur_35k_significant = sum(abs(eur_35k_z) > z_score_threshold),
      eur_31k_significant = sum(abs(eur_31k_z) > z_score_threshold),

      eur_35k_significant_p = eur_35k_significant/eur_31k_significant,
      multi_significant_p = multi_significant/eur_31k_significant,

      average_increase_eur_35k_compared_to_eur_31k = mean(eur_35k_z / eur_31k_z),
      average_increase_multi_compared_to_eur_31k = mean(multi_z / eur_31k_z))

  summary_statistics_diff %>% pivot_longer(cols=everything())

  # We can quantify: the number of primary eQTLs found, the
  summary_statistics_diff_cis_trans <- summary_statistics_merged_35k %>%
    ungroup() %>%
    group_by(typeCisTrans) %>%
    summarise(
      n = n(),
      all_singificant = sum(abs(z_score) > z_score_threshold),
      multi_significant = sum(abs(multi_z) > z_score_threshold),
      eur_35k_significant = sum(abs(eur_35k_z) > z_score_threshold),
      eur_31k_significant = sum(abs(eur_31k_z) > z_score_threshold),

      eur_35k_significant_p = eur_35k_significant/eur_31k_significant,
      multi_significant_p = multi_significant/eur_31k_significant,

      average_increase_eur_35k_compared_to_eur_31k = mean(eur_35k_z / eur_31k_z),
      average_increase_multi_compared_to_eur_31k = mean(multi_z / eur_31k_z))

  summary_statistics_diff_cis_trans %>% pivot_longer(cols=-typeCisTrans)


  meta_analysis_gait2 <- per_cohort_summary_statistics %>%
    filter(!is.na(beta)) %>%
    group_by(phenotype, variant_index) %>%
    mutate(weights = 1 / standard_error^2) %>%
    summarise(
      meta_beta = sum(beta * weights) / sum(weights),
      meta_se = sqrt(1 / sum(weights)),
      meta_z = meta_beta / meta_se,
      lco_beta = sum(beta[cohort != "GAIT2"] * weights[cohort != "GAIT2"]) / sum(weights[cohort != "GAIT2"]),
      lco_se = sqrt(1 / sum(weights[cohort != "GAIT2"])),
      lco_z = lco_beta / lco_se,
    )


  meta_analysis_gtex <- per_cohort_summary_statistics %>%
    filter(!is.na(beta), cohort != "GTEx_v9_EUR") %>%
    group_by(phenotype, variant_index) %>%
    mutate(weights = 1 / standard_error^2) %>%
    summarise(
      beta = sum(beta * weights) / sum(weights),
      se = sqrt(1 / sum(weights)))

  effect_lco <- meta_analysis_gait2 %>%
    inner_join(independent_variant_set_annot_filtered, by = c("variant_index", "phenotype")) %>%
    ungroup() %>%
    group_by(typeCisTrans) %>%
    summarise(
      all_more_significant = sum(abs(lco_z) < abs(meta_z)),
      all_more_significant_p = all_more_significant / n(),
      correlation = cor(lco_z, meta_z),
      ols = list(tidy(lm(meta_z ~ lco_z)))
    )

  ggplot(meta_analysis_gait2 %>%
           inner_join(independent_variant_set_annot_filtered, by = c("variant_index", "phenotype")), aes(lco_beta, meta_beta)) +
    geom_point(shape=16, size=0.5, alpha=0.5) +
    geom_smooth(method="lm")+
    geom_abline(color="salmon") +
    xlab("Betas (Meta-analysis w/o GAIT2)") +
    ylab("Betas (Meta-analysis with GAIT2)") +
    facet_grid(cols=vars(typeCisTrans)) + theme(aspect.ratio = 1)

  ggsave("leave_gait2_out_betas.pdf", width=10, height=5)



  ggplot(meta_analysis_gtex %>%
           inner_join(independent_variant_set_annot_filtered, by = c("variant_index", "phenotype")), aes(lco_beta, meta_beta)) +
    geom_point(shape=16, size=0.5, alpha=0.5) +
    geom_smooth(method="lm")+
    geom_abline(color="salmon") +
    xlab("Betas (Meta-analysis w/o GTEx)") +
    ylab("Betas (Meta-analysis with GTEx)") +
    facet_grid(cols=vars(typeCisTrans)) + theme(aspect.ratio = 1)

  ggsave("leave_gtex_out_betas.pdf", width=10, height=5)

  ggplot(meta_analysis_gtex %>%
           inner_join(independent_variant_set_annot_filtered, by = c("variant_index", "phenotype")), aes(lco_z, meta_z)) +
    geom_point(shape=16, size=0.5, alpha=0.5) +
    geom_smooth(method="lm")+
    geom_abline(color="salmon") +
    xlab("t-scores (Meta-analysis w/o GTEx)") +
    ylab("t-scores (Meta-analysis with GTEx)") +
    facet_grid(cols=vars(typeCisTrans)) + theme(aspect.ratio = 1)

  ggsave("leave_gtex_out_z.pdf", width=10, height=5)


  effect_lco_gtex <- meta_analysis_gtex %>%
    inner_join(independent_variant_set_annot_filtered, by = c("variant_index", "phenotype")) %>%
    ungroup() %>%
    group_by(typeCisTrans) %>%
    summarise(
      all_more_significant = sum(abs(lco_z) < abs(meta_z)),
      all_more_significant_p = all_more_significant / n(),
      correlation = cor(lco_z, meta_z),
      ols = list(tidy(lm(meta_z ~ lco_z)))
    )
  # Compared to the subset of 31K European-only samples,
  # adding all non-European samples increases sample size to 35988,
  # increases z-scores 1.05 times, and the number of significant primary eQTLs 1.08 times.

  # Compared to the subset of 31K European-only samples,
  # adding European samples to 36017 increases z-scores 1.15 times and the number of significant primary eQTLs 1.30 times
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}