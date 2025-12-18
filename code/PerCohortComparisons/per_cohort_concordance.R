#!/usr/bin/env Rscript
#' Generate plots for replication results

# Load libraries
library(tidyverse)
library(arrow)
library(data.table)

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

pop_colors <- c("African (Tanzania)" = "red",
                "African (African American)" = "tomato",
                "European (Ashkenazi)" = "darkturquoise",
                "European (Finnish)" = "steelblue4",
                "European (Estonian)" = "steelblue1",
                "European (West)" = "lightblue",
                "European (West and South)" = "lightblue3",
                "East Asian (Japanese)" = "forestgreen",
                "African (Moroccan)" = "deeppink3",
                "Middle-Eastern (Qatari)" = "pink4",
                "European (Toscani)" = "lightblue4",
                "European (Spanish)" = "deepskyblue3",
                "South Asian (Bangladeshi)" = "orange",
                "South Asian (Indian)" = "tan2",
                "East Asian (Malaysia)" = "darkolivegreen2",
                "East Asian (Chinese)" = "springgreen3")


# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))
  location <- str_split(variables["available_genes"], "freeze3")[[1]][1]

  cohort_info <- fread("") # Here goes the per cohort information
  master_table <- fread(variables["master_table"]) %>% select(cohort, cohort_new_name, dataset_n, dataset_platform, dataset_ancestry)

  # Process input
  variant_reference <- read_parquet(variables["variant_reference"])

  # Load the main file of independent variants
  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type)) %>%
    select(variant_index, phenotype, beta, standard_error, type, typeCisTrans) %>%
    mutate(z_score = beta / standard_error)

  per_cohort_extractions <- "" # Here go the per-cohort results
  per_cohort_summary_statistics <- open_dataset(per_cohort_extractions) %>% collect()
  per_cohort_summary_statistics <- per_cohort_summary_statistics %>%
    mutate(z_score = beta / standard_error)

  summary_statistics_merged <- independent_variant_set_annot_filtered %>%
    left_join(per_cohort_summary_statistics, by = c("variant_index", "phenotype"), suffix=c("_meta", "_cohort")) %>%
    mutate(effect_concordant = if_else(sign(z_score_meta) == sign(z_score_cohort), "concordant", "discordant")) %>%
    left_join(master_table %>% select(cohort, cohort_new_name), by = c("cohort" = "cohort"))

  statistics_summarised <- summary_statistics_merged %>%
    group_by(cohort, cohort_new_name, type) %>%
    summarise(
      n_total = n(),
      n_cohort = sum(!is.na(z_score_cohort)),
      n_concordant = sum(effect_concordant == "concordant", na.rm = TRUE),
      pct_plotted = 100 * n_cohort / n_total,
      pct_concordant = 100 * n_concordant / n_cohort,
      lm_summary = nest(tidy(lm(beta_cohort ~ beta_meta))),
    ) %>%
    unnest(cols = c(lm_summary)) %>%
    unnest_wider(col=data) %>%
    unnest_wider(estimate, names_sep="_") %>%
    mutate(
      label = sprintf(
        "N plotted: %d (%.1f%%)\nConcordant: %d (%.1f%%)\nSlope of betas: %.2f",
        n_total, pct_plotted, n_concordant, pct_concordant, estimate_2
      )
    )

  library(ggplot2)

  # Define colors
  vermillion <- "#D55E00"
  grey_col <- "grey40"

  # Main plot
  p <- ggplot(summary_statistics_merged %>% filter(type == "cis"), aes(x = z_score_meta, y = z_score_cohort)) +
    geom_hline(yintercept = 0, color="grey90", linewidth=(0.5 / (ggplot2::.pt * 72.27 / 96))) +
    geom_vline(xintercept = 0, color="grey90", linewidth=(0.5 / (ggplot2::.pt * 72.27 / 96))) +
    rasterize(geom_point(
      aes(color = effect_concordant),
      alpha = 0.2, size = 0.2, shape=16
    ), dpi=300) +
    scale_color_manual(
      values = c("concordant" = grey_col, "discordant" = vermillion),
      name = "Effect direction"
    ) +
    geom_text(
      data = statistics_summarised %>% filter(type == "cis"),
      aes(x = Inf, y= -Inf, label=label),
      inherit.aes=F,
      hjust="inward",
      vjust="inward",
      size = 4 / .pt, color = "black"
    ) +
    facet_wrap(
      ~ cohort_new_name,
      nrow = 9,
      scales = "free",   # same axis scales for all
      labeller = label_wrap_gen(width = 15)  # multi-line facet labels
    ) +
    labs(
      x = "Meta-analysis Z-score",
      y = "Per-cohort Z-score",
      title = "Comparison of meta-analyzed and per-cohort Z-scores of cis-eQTLs"
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.clip = "off",
      strip.text = element_text(hjust=0)
    )


  ggsave("per_cohort_concordance_cis.pdf", p, width=180, height=267, units="mm")


  # Main plot
  p <- ggplot(summary_statistics_merged %>% filter(type == "trans"), aes(x = z_score_meta, y = z_score_cohort)) +
    geom_hline(yintercept = 0, color="grey90", linewidth=(0.5 / (ggplot2::.pt * 72.27 / 96))) +
    geom_vline(xintercept = 0, color="grey90", linewidth=(0.5 / (ggplot2::.pt * 72.27 / 96))) +
    rasterize(geom_point(
      aes(color = effect_concordant),
      alpha = 0.2, size = 0.2, shape=16
    ), dpi=300) +
    scale_color_manual(
      values = c("concordant" = grey_col, "discordant" = vermillion),
      name = "Effect direction"
    ) +
    geom_text(
      data = statistics_summarised %>% filter(type == "trans"),
      aes(x = Inf, y= -Inf, label=label),
      inherit.aes=F,
      hjust="inward",
      vjust="inward",
      size = 4 / .pt, color = "black",
    ) +
    facet_wrap(
      ~ cohort_new_name,
      nrow = 9,
      scales = "free",   # same axis scales for all
      labeller = label_wrap_gen(width = 15)  # multi-line facet labels
    ) +
    labs(
      x = "Meta-analysis Z-score",
      y = "Per-cohort Z-score",
      title = "Comparison of meta-analyzed and per-cohort Z-scores of trans-eQTLs"
    ) +
    theme(
      legend.position = "bottom",
      strip.clip = "off",
      strip.text = element_text(hjust=0)
    )


  ggsave("per_cohort_concordance_trans.pdf", p, width=180, height=267, units="mm")



}

if (sys.nframe() == 0 && !interactive()) {
  main()
}