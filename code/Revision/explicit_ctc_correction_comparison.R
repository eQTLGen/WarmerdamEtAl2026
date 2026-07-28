#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)
library(broom)
library(bit64)


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
  # Get resources
  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))

  variant_classes <- fread(variables["variant_coloc_classes"]) %>% select(phenotype, variant_index, variant_type) %>%
    mutate(coloc = variant_type != "Not GWAS colocalising")

  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type)) %>%
    select(variant_index, phenotype, typeCisTrans)

  # Process input
  eqtls <- fread(variables["independent_variant_set"]) %>%
    mutate(variant_index = as.integer64(variant_index))

  ds_corrected <- open_dataset("eQTLGen/cell_type_composition/eqtl_mapping/output/cell_type_corrected_2026-06-09/eqtls/meta")
  ds_uncorrected <- open_dataset("eQTLGen/cell_type_composition/eqtl_mapping/output/cell_type_uncorrected_2026-06-09/eqtls/meta")
  ds_corrected_noExpPCs <- open_dataset("eQTLGen/cell_type_composition/eqtl_mapping/output/cell_type_corrected_noExpPCs_2026-06-11/eqtls/meta")

  eqtls_adj <- ds_corrected %>% semi_join(eqtls, by = c("variant_index", "phenotype")) %>% collect()
  eqtls_unadj <- ds_uncorrected %>% semi_join(eqtls, by = c("variant_index", "phenotype")) %>% collect()
  eqtls_adj_noExpPCs <- ds_corrected_noExpPCs %>% semi_join(eqtls, by = c("variant_index", "phenotype")) %>% collect()

  saveRDS(eqtls_adj_noExpPCs, "eQTLs_explicitly_adjusted_for_ctc_noExpPCs_lead_variants_20260522_20260611.rds")
  saveRDS(eqtls_adj, "eQTLs_explicitly_adjusted_for_ctc_lead_variants_20260522_20260611.rds")
  saveRDS(eqtls_unadj, "eQTLs_not_explicitly_adjusted_for_ctc_lead_variants_20260522_20260609.rds")

  eqtls_adj_noExpPCs <- readRDS( "eQTLs_explicitly_adjusted_for_ctc_noExpPCs_lead_variants_20260522_20260611.rds")
  eqtls_adj <- readRDS( "eQTLs_explicitly_adjusted_for_ctc_lead_variants_20260522_20260611.rds")
  eqtls_unadj <- readRDS("eQTLs_not_explicitly_adjusted_for_ctc_lead_variants_20260522_20260609.rds")

  joined <- inner_join(eqtls_adj, eqtls_unadj, by =c("variant_index", "phenotype", "sample_size"), suffix=c("_adj", "_unadj")) %>%
    mutate(t_adj = beta_adj / standard_error_adj,
           t_unadj = beta_unadj / standard_error_unadj)  %>%
    left_join(variant_classes, by = c("variant_index", "phenotype")) %>%
    inner_join(independent_variant_set_annot_filtered, by = c("variant_index", "phenotype"))

  joined %>% summarise(
    adj_larger = sum(beta_adj > beta_unadj),
    unadj_larger = sum(beta_unadj > beta_adj)
  )

  correlations <- joined %>% group_by(typeCisTrans) %>%
    summarise(tidy(cor.test(beta_adj, beta_unadj)))

  joined %>% filter(typeCisTrans == "trans") %>%
    reframe(lm_res = tidy(lm(beta_adj ~ beta_unadj))) %>%
    unnest(lm_res) %>%
    filter(term == "beta_unadj") %>%
    mutate(
      t = (estimate - 1) / std.error,
      conf.upper = estimate + (qnorm(0.975, lower.tail = T) * std.error),
      conf.lower = estimate - (qnorm(0.975, lower.tail = T) * std.error),
      res = sprintf("%.3f (%.3f-%.3f)", estimate, conf.lower, conf.upper))

  joined %>% filter(typeCisTrans == "trans") %>%
    group_by(variant_type) %>%
    reframe(lm_res = tidy(lm(beta_adj ~ beta_unadj))) %>%
    unnest(lm_res) %>%
    filter(term == "beta_unadj") %>%
    mutate(t = (estimate - 1) / std.error, res = sprintf("%.6f", estimate))

  ggplot(joined, aes(beta_unadj, beta_adj, color = variant_type)) +
    geom_abline(color = "black") +
    geom_point(size=0.5, shape=16, alpha=0.5) +
    geom_smooth(method="lm", fullrange=T, linewidth=(0.5 / (ggplot2::.pt * 72.27 / 96)), se=F) +
    ylab("betas, analysis 2\n(20 expression PCs + direct adjustment for CTC measurements)") +
    xlab("betas, analysis 1 (20 expression PCs only)") +
    theme(legend.position = "bottom") +
    facet_wrap(vars(typeCisTrans), scales="free")

  ggsave("comparison_joined_20PCsOnly_vs_20PCsAndCtc.pdf", width=168, height=96, units="mm")

  ggplot(joined, aes(t_unadj, t_adj, color = variant_type)) +
    geom_abline(color = "black") +
    geom_point(size=0.5, shape=16, alpha=0.8) +
    geom_smooth(method="lm", fullrange=T) +
    facet_wrap(vars(typeCisTrans), scales="free")
  # Work



  joined_2 <- inner_join(eqtls_unadj, eqtls_adj_noExpPCs, by =c("variant_index", "phenotype", "sample_size"), suffix=c("_unadj", "_adj_noExpPCs")) %>%
    mutate(t_unadj = beta_unadj / standard_error_unadj,
           t_adj_noExpPCs = beta_adj_noExpPCs / standard_error_adj_noExpPCs) %>%
    left_join(variant_classes, by = c("variant_index", "phenotype")) %>%
    inner_join(independent_variant_set_annot_filtered, by = c("variant_index", "phenotype"))

  joined_2 %>% filter(typeCisTrans == "trans") %>%
    group_by(variant_type) %>%
    reframe(lm_res = tidy(lm(beta_adj_noExpPCs ~ beta_unadj))) %>%
    unnest(lm_res) %>%
    filter(term == "beta_unadj") %>%
    mutate(
      t = (estimate - 1) / std.error,
      conf.upper = estimate + (qnorm(0.975, lower.tail = T) * std.error),
      conf.lower = estimate - (qnorm(0.975, lower.tail = T) * std.error),
      res = sprintf("%.3f (%.3f-%.3f)", estimate, conf.lower, conf.upper))

  joined_2 %>% group_by(typeCisTrans) %>%
    summarise(tidy(cor.test(beta_adj_noExpPCs, beta_unadj)))

  ggplot(joined_2, aes(beta_unadj, beta_adj_noExpPCs, color = variant_type)) +
    geom_abline(color = "black") +
    geom_point(size=0.5, shape=16, alpha=0.5) +
    geom_smooth(method="lm", fullrange=T, linewidth=(0.5 / (ggplot2::.pt * 72.27 / 96)), se=F) +
    ylab("betas, analysis 3\n(no expression PCs, direct adjustment for CTC measurements)") +
    xlab("betas, analysis 1 (20 expression PCs only)") +
    theme(legend.position = "bottom") +
    facet_wrap(vars(typeCisTrans), scales="free")

  ggsave("comparison_joined_20PCsOnly_vs_CtcOnly.pdf", width=168, height=96, units="mm")

  ggplot(joined, aes(t_unadj, t_adj_noExpPCs, color = variant_type)) +
    geom_abline(color = "black") +
    geom_point(size=0.5, shape=16, alpha=0.8) +
    geom_smooth(method="lm", fullrange=T, linewidth=(0.5 / (ggplot2::.pt * 72.27 / 96))) +
    facet_wrap(vars(typeCisTrans), scales="free")
  # Work
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}