#!/usr/bin/env Rscript
#' Generate plots for replication results

# Load libraries
library(tidyverse)
library(arrow)
library(data.table)
library(ggrepel)
library(patchwork)
library(broom)
library(ggrastr)

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

#' Maximum Independent Set of Uncorrelated Genes
#'
#' This function identifies the maximum set of uncorrelated genes (independent vertices)
#' from a boolean adjacency matrix. It uses an algorithm that selects the vertex with the
#' least amount of edges and removes all vertices that are connected to it.
#'
#' @param adjacency_matrix A boolean adjacency matrix representing gene correlations (1 if correlated, 0 if not).
#' @param nodes A character vector containing the names of the genes corresponding to the rows/columns of the matrix.
#' @return A character vector of genes that are part of the maximum independent set.
#' @examples
#' adj_matrix <- matrix(c(1, 1, 0, 1, 1, 1, 0, 0, 1), nrow = 3, byrow = TRUE)
#' gene_names <- c("Gene1", "Gene2", "Gene3")
#' maximum_independent_set(adj_matrix, gene_names)
maximum_independent_set <- function(adjacency_matrix, nodes) {
  # Initialize the set of independent vertices to be empty
  independent_vertices <- character(0)
  active_nodes <- nodes == nodes
  names(active_nodes) <- nodes

  colnames(adjacency_matrix) <- nodes
  rownames(adjacency_matrix) <- nodes

  diag(adjacency_matrix) <- T

  # First collect every gene that is totally independent.
  # Removing these will not have an effect on other genes.
  unconnected_vertices <- rowSums(adjacency_matrix) == 1
  independent_vertices <- c(independent_vertices, nodes[unconnected_vertices])

  with_edges <- colSums(adjacency_matrix[unconnected_vertices, active_nodes] == T) > 0
  active_nodes[with_edges] <- F

  print(length(unconnected_vertices))
  print(sum(active_nodes))

  # Loop while the matrix is not empty
  while (any(active_nodes)) {
    # Find the vertex with the least amount of edges
    least_connected <- names(base::which.min(colSums(adjacency_matrix[active_nodes, active_nodes, drop=F])))
    cat(least_connected, "\n")

    # Add the least connected vertex to the independent set
    independent_vertices <- c(independent_vertices, least_connected)

    # Remove all vertices that have an edge with the least connected vertex
    with_edges <- adjacency_matrix[least_connected,][active_nodes, drop=F] == T

    active_nodes[names(with_edges)[with_edges]] <- F

    cat(sum(active_nodes), "\n")
  }

  return(independent_vertices)
}


# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  type <- "cis"

  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))
  location <- str_split(variables["available_genes"], "freeze3")[[1]][1]

  cohort_info <- fread("cohort_info.txt")
  master_table <- fread(variables["master_table"]) %>% select(cohort, cohort_new_name, dataset_n, dataset_platform, dataset_ancestry)

  # Process input
  variant_reference <- read_parquet(variables["variant_reference"])

  # Load the main file of independent variants
  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type)) %>%
    select(variant_index, phenotype, beta, standard_error, type, typeCisTrans, SusieRss_CS) %>%
    mutate(z_score = beta / standard_error)

  per_cohort_extractions <- variables["per_cohort_extractions"]
  per_cohort_summary_statistics <- readRDS(per_cohort_extractions)
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
      cor_summary = nest(tidy(cor.test(z_score_meta, z_score_cohort)))
    ) %>%
    unnest(cols = c(cor_summary)) %>%
    unnest_wider(col=data) %>%
    mutate(
      label = sprintf(
        "N plotted: %d (%.1f%%)\nConcordant: %d (%.1f%%)\nCorrelation: %.2f",
        n_total, pct_plotted, n_concordant, pct_concordant, estimate
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


  ggsave("per_cohort_concordance_cis_20260623.pdf", p, width=180, height=267, units="mm")


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


  ggsave("per_cohort_concordance_trans_20260623.pdf", p, width=180, height=267, units="mm")

  variants <- fread(variables["uncorrelated_variants"], header=F, col.names = c("variants"))$variants

  variants_uncorrelated <- variant_reference %>% filter(variant %in% variants)

  uncorrelated_genes <- readRDS(variables["uncorrelated_genes"])

  ld_matrices <- readRDS(variables["finemapped_variant_ld"])

  inclusion_table <- read_parquet(variables["variant_inclusion_table"])
  inclusion_mat <- inclusion_table %>% as.data.table() %>% as.matrix(rownames ="variant_index")
  variants_available <- rownames(inclusion_mat)[rowSums(inclusion_mat) == 52]

  uncorrelated_fine_mapped_snps <- mapply(function(ld_matrix) {
    str(ld_matrix)
    variants_available_this_chr <- variants_available[variants_available %in% rownames(ld_matrix)]
    ld_matrix_complete_vars <- ld_matrix[variants_available_this_chr, variants_available_this_chr]
    str(ld_matrix_complete_vars)
    maximum_independent_set(ld_matrix_complete_vars^2 > 0.01, rownames(ld_matrix_complete_vars))
  }, ld_matrices, SIMPLIFY=F)

  uncorrelated_fine_mapped_snps_bound <- as.numeric(unlist(uncorrelated_fine_mapped_snps))

  per_cohort_summary_statistics_uncorrelated <- per_cohort_summary_statistics %>%
    filter(phenotype %in% uncorrelated_genes, variant_index %in% uncorrelated_fine_mapped_snps_bound)

  betas_uncorrelated_variants_genes <- per_cohort_summary_statistics_uncorrelated %>%
    pivot_wider(id_cols = c(variant_index, phenotype), names_from = "cohort", values_from = "beta") %>%
    group_by(variant_index) %>%
    slice_sample(n=1) %>%
    group_by(phenotype) %>%
    slice_sample(n=1) %>%
    ungroup()

  mat_betas <- betas_uncorrelated_variants_genes %>%
    select(-variant_index, -phenotype) %>%
    as.data.table() %>% as.matrix()

  sum(is.na(mat_betas)) / length(mat_betas) * 100

  mat_betas_filled <- t(apply(mat_betas, 1, function(row) {
    row[is.na(row)] <- mean(row, na.rm = TRUE)
    return(row)
  }))

  mat_betas_scaled <- t(scale(mat_betas_filled, center = T, scale =T))
  svd.res <- svd(mat_betas_scaled)

  coords <- as.data.table(mat_betas_scaled %*% svd.res$v, keep.rownames = "cohort")
  colnames(coords)[2:5] <- c("PC1", "PC2", "PC3", "PC4")
  var_explained <- svd.res$d^2 / sum(svd.res$d^2)

  mat_betas_correlations <- cor(mat_betas_filled)

  as.data.table(mat_betas_correlations, keep.rownames = "Cohort1") %>%
    pivot_longer(cols=-Cohort1, values_to = "Cor", names_to = "Cohort2") %>%
    ggplot(aes(Cohort1, Cohort2, fill=Cor)) +
    geom_tile() + scale_fill_gradient2(low="darkblue", high="darkgreen", guide="colorbar", mid = 0)

  ggsave("cohort-cohort_cor.pdf")

  eig <- eigen(mat_betas_correlations)

  coords <- as.data.frame(eig$vectors %*% diag(sqrt(eig$values)))
  colnames(coords)[1:ncol(coords)] <- paste0("PC", 1:ncol(coords))

  coords$cohort <- rownames(mat_betas_correlations)
  var_explained <- eig$values / sum(eig$values)

  coords_annotated <- coords %>% left_join(master_table, by = c("cohort")) %>% left_join(cohort_info, by = c("cohort_new_name" = "cohort"))

  p1 <- ggplot(coords_annotated, aes(x = -PC1, y = PC2, label = cohort)) +
    geom_point(aes(fill = ancestry, shape = exp_platform, size=dataset_n)) +
    geom_text_repel(data = . %>% filter(cohort == "CAD"), size = 2) +
    scale_fill_manual(values = pop_colors) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    labs(
      title = "PCA of pruned per-cohort results",
      x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
    ) + guides(fill = guide_legend(title="Ancestry", override.aes=list(shape=21)), size = guide_legend(title="Sample size"), shape = guide_legend(title="Expression platform")) + theme(aspect.ratio = 1)

  p2 <- ggplot(coords_annotated, aes(x = -PC3, y = -PC4, label = cohort)) +
    geom_point(aes(fill = ancestry, shape = exp_platform, size=dataset_n)) +
    geom_text_repel(data = . %>% filter(cohort == "CAD"), size = 2) +
    scale_fill_manual(values = pop_colors) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    labs(
      x = paste0("PC3 (", round(var_explained[3] * 100, 1), "%)"),
      y = paste0("PC4 (", round(var_explained[4] * 100, 1), "%)")
    ) + guides(fill = guide_legend(title="Ancestry", override.aes=list(shape=21)), size = guide_legend(title="Sample size"), shape = guide_legend(title="Expression platform")) + theme(aspect.ratio = 1)

  p_patch <- p1 + p2 + plot_layout(ncol = 1, guides = 'collect')

  ggsave("pca_coords.pdf", p_patch, height=8, width=6)


  # Association tests:

  coords_annotated %>%
    group_modify(~ {tidy(lm(.x$dataset_n ~ .x$exp_platform + .x$PC1 + .x$PC2))})

  coords_annotated %>%
    group_by(exp_platform) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    group_modify(~ {tidy(oneway.test(.x$PC2 ~ .x$exp_platform, var.equal=F))})

  coords_annotated %>%
    mutate(PC1_no_exp = residuals(lm(PC1 ~ exp_platform))) %>%
    ungroup() %>%
    group_modify(~ {tidy(cor.test(.x$PC1_no_exp, .x$dataset_n))})

  mds_coordinates <- fread(variables["mds_coordinates"]) %>%
    mutate(Cohort = cohort_name_mapping[Cohort])

  coords_annotated_ancestry <- coords_annotated %>%
    inner_join(mds_coordinates, by = c("cohort" = "Cohort"))

  library(psych)

  pairwise_correlations <- corr.test(coords_annotated_ancestry[, paste0("PC", 1:4)], coords_annotated_ancestry[, paste0("MDS", 1:4)])

  lm.res <- coords_annotated_ancestry %>%
    group_modify(~ {tidy(lm(as.matrix(.x[paste0("PC", 1:4)]) ~ MDS1 + MDS2 + MDS3 + MDS4 + exp_platform + geno_platform + dataset_n, data=.x))})


  res_anova <- lapply(paste0("PC", 1:4), function(pc) {
    fit <- lm(get(pc) ~ MDS1 + MDS2 + MDS3 + MDS4 + exp_platform + geno_platform + dataset_n, data = coords_annotated_ancestry)
    a <- anova(fit)

    ss_total <- sum(a$"Sum Sq")
    a$prop_var <- a$"Sum Sq" / ss_total
    a
    return(a)
  })

  names(res_anova) <- paste0("PC", 1:4)
  anova_res_tib <- bind_rows(res_anova, .id = "PC")

  anova_res_tib %>%
    filter(str_starts(term, "MDS")) %>%
    group_by(PC) %>%
    summarise(sum(prop_var))

  res_r2 <- lapply(paste0("PC", 1:4), function(pc) {

    full <- lm(get(pc) ~ MDS1 + MDS2 + MDS3 + MDS4 +
      exp_platform + geno_platform + dataset_n,
               data = coords_annotated_ancestry)

    r2_full <- summary(full)$r.squared

    vars <- attr(terms(full), "term.labels")

    sapply(vars, function(v) {
      reduced <- update(full, paste(". ~ . -", v))
      r2_full - summary(reduced)$r.squared
    })
  })

  res_r2 <- do.call(rbind, res_r2)


  per_cohort_summary_statistics_proc <- per_cohort_summary_statistics %>%
    filter(!is.na(z_score)) %>%
    inner_join(independent_variant_set_annot_filtered %>%
                 filter(!is.na(SusieRss_CS)) %>%
                 select(typeCisTrans, variant_index, phenotype), by=c("phenotype", "variant_index")) %>%
    mutate(p = 2*pnorm(abs(z_score), lower.tail = F),
           significant = p<0.05)

  single_cohort_significant <- per_cohort_summary_statistics_proc %>%
    group_by(phenotype, variant_index) %>%
    summarise(n_nominally_significant_cohorts = sum(significant)) %>%
    filter(n_nominally_significant_cohorts == 1)

  single_cohort_significant

  single_cohort_significant_after_lco <- per_cohort_summary_statistics_proc %>%
    inner_join(single_cohort_significant) %>%
    group_by(phenotype, variant_index, typeCisTrans) %>%
    mutate(weights = 1 / standard_error^2) %>%
    summarise(
      sigc_beta = beta[significant],
      sigc_se = standard_error[significant],
      significant_cohort = cohort[significant],
      meta_beta = sum(beta * weights) / sum(weights),
      meta_se = sqrt(1 / sum(weights)),
      meta_z = meta_beta / meta_se,
      meta_p = 2*pnorm(abs(meta_z), lower.tail = F),
      meta_beta_lco = sum(beta[!significant] * weights[!significant]) / sum(weights[!significant]),
      meta_se_lco = sqrt(1 / sum(weights[!significant])),
      meta_z_lco = meta_beta_lco / meta_se_lco,
      meta_p_lco = 2*pnorm(abs(meta_z_lco), lower.tail = F),
      meta_sig = meta_p < 5e-8, meta_sig_lco = meta_p_lco < 5e-8,

      z_diff = (sigc_beta - meta_beta_lco) /
        sqrt(sigc_se^2 + meta_se_lco^2),
      p_diff = 2 * pnorm(abs(z_diff), lower.tail = FALSE),
    ) %>% filter(meta_sig)

  single_cohort_significant_after_lco %>% group_by(significant_cohort, typeCisTrans) %>% count()

  single_cohort_significant_after_lco %>% select(meta_p, meta_p_lco, meta_sig, meta_sig_lco, typeCisTrans) %>%
    filter(meta_sig) %>%
    group_by(meta_sig, meta_sig_lco, typeCisTrans, meta_p_lco < 0.05) %>%
    tally()
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}