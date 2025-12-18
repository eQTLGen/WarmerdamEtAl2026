#!/usr/bin/env Rscript
#' This script calculates the results of the manuscript up to the replication analyses
#' - Number of cis- and trans-eQTLs
#' - Number of Credible sets
#' - Relationship between the number if credible sets, gene expression and gene conservation
#' - Variant class enrichments
#' - VEP enrichments


# Load libraries
library(data.table)
library(arrow)
library(tidyverse)
library(dtplyr)
library(dbplyr)
library(scales)
library(ggh4x)
library(rtracklayer)
library(ggpubr)
library(patchwork)
library(RobustRankAggreg)
library(broom)
library(fgsea)

# Declare constants
old <- theme_set(theme_classic(base_size = 8))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  panel.grid.minor = element_blank(),
  text = element_text(family="Helvetica"),
  axis.title = element_text(colour="black"),
  title = element_text(colour = "#595A5C", face = "plain")
)

eqtl_overview_plot <- function(data) {

  # Define color ramp: light to base color
  base_color <- "#0072B2"
  # A strong blue
  light_color <- "#C3E9FD"
  # A lighter blue (or use "white")

  # Create function for n shades from light to base
  n_colors <- 11
  # Adjust based on your L range
  color_func <- colorRampPalette(c(base_color, light_color))

  types <- c("cis", "trans")
  plot_pairs <- list()

  for (t in types) {
    data_sub <- filter(data, typeCisTrans == t)

    # Bar plot
    bar_plot <- ggplot(data_sub, aes(bin, fill = n_CSs_per_gene_ord_capped)) +
      geom_bar(position = "fill") +
      scale_fill_manual(name = "Number of eQTLs", values = c("grey90", color_func(11)), guide=NULL) +
      scale_x_discrete(name = "Gene expression in blood",
                       breaks = c("Bin1", "Bin10"),
                       labels = c("Low", "High")) +
      scale_y_continuous(name="Proportion of genes\nshowing eQTL effect")
      theme(legend.position = "none")

    # Correlation data
    cor_data <- data_sub %>%
      group_by(bin) %>%
      nest() %>%
      summarise(
        Cor = map(data, ~ tidy(cor.test(.x$n_CSs_per_gene, .x$loeuf_mane, method = "spearman"))),
        MWU = map(data, ~ tidy(wilcox.test(.x$pLI_mane[.x$any_CS], .x$pLI_mane[!.x$any_CS])))
      ) %>%
      unnest(c("Cor", "MWU"), names_sep=":") %>%
      pivot_longer(cols=ends_with("p.value"), names_sep=":", names_to=c("test", "statistic")) %>%
      mutate(
        p_lab = if_else(value > 0.05, "NS", gsub("e", " %*% 10^", sprintf("%.0e", value))))

    stats <- ggplot(cor_data, aes(x=1, y=test, label=p_lab)) +
      geom_text(size=4, size.unit = "pt", parse=T) + theme_minimal() +
      scale_y_discrete(breaks=c("MWU", "Spearman")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
           panel.border = element_blank(), axis.text.x = element_blank(), strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.text.y = element_text(angle=45, size=1),
            plot.margin = unit(c(0, 0, 0, 0), "null"),
           axis.ticks = element_blank()) + xlab(NULL) + ylab(NULL) +
      coord_cartesian(clip="off") +
      facet_wrap(vars(bin), nrow = 1)

      # Correlation plot
    cor_plot <- ggplot(data_sub, aes(n_CSs_per_gene, loeuf_mane, color = n_CSs_per_gene_ord_capped)) +
      stat_summary(fun.data = mean_se, fun.args = list(mult = 1), shape = 20,
                   size=(0.25 / (ggplot2::.pt * 72.27/96)), linewidth=(0.5 / (ggplot2::.pt * 72.27/96)),
                   geom = "pointrange") +
      scale_color_manual(name = "Number of eQTLs", values = c("grey90", color_func(11))) +
      scale_x_continuous(name = "number of eQTLs", limits = c(0, 11), breaks = c(0, 10)) +
      coord_cartesian(clip = "on") +
      scale_y_continuous(name = "LOEUF", breaks = c(0.5, 0.75, 1.0, 1.25)) +
      geom_smooth(method = "lm", aes(color = -1), color = "#D55E00",
                  linewidth=(1 / (ggplot2::.pt * 72.27/96))) +
      facet_wrap(vars(bin), nrow = 1) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.text.x= element_text(size=4),
            legend.position = "none")

    violins <- ggplot(data_sub,
           aes(x=any_CS, y=pLI_mane, fill=any_CS)) +
      geom_violin(linewidth=0, color="white") +
      stat_summary(fun = "median", geom = "crossbar", color = "#D55E00", linewidth=(0.5 / (ggplot2::.pt * 72.27/96))) +
      scale_fill_manual(values=c("grey90", color_func(11)[6])) +
      scale_y_continuous(name="pLI score", breaks=c(0, 1)) +
      scale_x_discrete(name="eGene", breaks=c("FALSE", "TRUE"), labels=c("no", "yes")) +
      coord_cartesian(ylim = c(0, 1.0), clip = "off") +
      facet_grid(cols=vars(typeCisTrans)) +
      facet_wrap(vars(bin), nrow = 1) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.text.x= element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x= element_blank(),
            legend.position = "none")

      plot_pairs[[t]] <- list(violins, stats, cor_plot, bar_plot)
  }

  layout <- "
  AAAAAAAABBBBBBBB
  CCCCCCCCDDDDDDDD
  EEEEEEEEFFFFFFFF
  GGGGGGGGHHHHHHHH
  #IIIII###JJJJJ##
  "

  final_plot <- plot_pairs$cis[[1]] + plot_pairs$trans[[1]] +
    plot_pairs$cis[[2]] + plot_pairs$trans[[2]] +
    plot_pairs$cis[[3]] + plot_pairs$trans[[3]] +
    plot_spacer() + plot_spacer() +
    plot_pairs$cis[[4]] + plot_pairs$trans[[4]] +
    plot_layout(heights = c(0.5, 0.2, 0.8, 0.1, 1.5), design=layout)

  return(final_plot)
}

eqtl_credible_set_leouf_supplement <- function(per_gene_type_combination_cis) {

  # Define color ramp: light to base color
  base_color <- "#0072B2"
  # A strong blue
  light_color <- "#C3E9FD"
  # A lighter blue (or use "white")

  # Create function for n shades from light to base
  n_colors <- 11
  # Adjust based on your L range
  color_func <- colorRampPalette(c(base_color, light_color))

  # Correlation data
  cor_data <- per_gene_type_combination_cis %>%
    group_by(bin) %>%
    nest() %>%
    summarise(
      Cor = map(data, ~ tidy(cor.test(.x$n_CSs_per_gene, .x$loeuf_mane, method = "spearman"))),
    ) %>%
    unnest(c("Cor")) %>%
    mutate(
      p_lab = if_else(p.value > 0.05, "NS", gsub("e", " %*% 10^", sprintf("%.0e", p.value))))

  res <- per_gene_type_combination_cis %>% ungroup() %>% mutate(n_CSs_per_gene_ord = ordered(n_CSs_per_gene))

  violin_dotplot <- ggplot(res,
                        aes(x=n_CSs_per_gene, y=loeuf_mane, fill = n_CSs_per_gene_ord_capped, group = n_CSs_per_gene_ord)) +
    geom_violin(linewidth = 0) +
    geom_boxplot(width=0.1, linewidth = 0.2, median.linewidth = 1, outlier.size = 0.2, fill="white", median.color = "black") +
    labs(x = "Number of eQTL credible sets") +
    scale_fill_manual(name = "Number of eQTLs", values = c("grey90", color_func(11))) +
    scale_x_continuous(name = "Number of eQTLs", limits = c(-0.5, 10.5), breaks = c(0, 10)) +
    scale_y_continuous(name = "LOEUF", breaks = c(0, 0.5, 1.0, 1.5, 2)) +
    geom_smooth(aes(x = n_CSs_per_gene, y = loeuf_mane), method = "lm", color = "#D55E00", inherit.aes=F,
                linewidth=(1 / (ggplot2::.pt * 72.27/96))) +
    coord_cartesian(clip=F) +
    geom_text(data = cor_data, aes(x=10.5, y=0.05, label=p_lab), color = "#D55E00",
              hjust = "inward", size=2,
              vjust = "inward", parse=T, inherit.aes=F, label.padding = unit(1, "lines"),) +
    facet_wrap(vars(bin), ncol = 4) +
    theme(legend.position = "bottom")

  return(violin_dotplot)
}

plot_beta_lof_summary_split_facets <- function(data, cis_trans_color_scale) {
  # Pre-calculate summaries
  beta_lof_z_summary <- data %>%
    filter(any_CS, n_CSs_per_gene > 0, !is.na(lof_z_mane), is.finite(lof_z_mane)) %>%
    group_by(typeCisTrans) %>%
    nest() %>%
    summarise(
      lof_beta_cor = map(data, ~ tidy(cor.test(.x$lof_z_mane, log10(abs(.x$beta))))),
      median_abs_beta = map(data, ~ median(abs(.x$beta)))
    ) %>%
    unnest(c(lof_beta_cor, median_abs_beta)) %>%
    mutate(
      log_cor_label = paste0(sprintf("rho==%.2f", estimate), "*", format_pvalues(p.value))
    )

  beta_summary <- data %>%
    filter(any_CS, n_CSs_per_gene > 0, !is.na(lof_z_mane), is.finite(lof_z_mane)) %>%
    nest() %>%
    summarise(
      MWU = map(data, ~ {
        df <- .x
        cis_betas <- abs(df$beta[df$typeCisTrans == "cis"])
        trans_betas <- abs(df$beta[df$typeCisTrans == "trans"])
        if (length(cis_betas) > 0 && length(trans_betas) > 0) {
          tidy(wilcox.test(cis_betas, trans_betas))
        } else {
          tibble(statistic = NA, p.value = NA)
        }
      })
    ) %>%
    unnest(MWU) %>%
    mutate(mwu_label = format_pvalues(p.value))

  # Prepare data subsets
  cis_data <- data %>% filter(any_CS, typeCisTrans == "cis")
  trans_data <- data %>% filter(any_CS, typeCisTrans == "trans")

  x_min <- min(abs(cis_data$beta), abs(trans_data$beta))
  x_max <- max(abs(cis_data$beta), abs(trans_data$beta))
  y_min <- min(cis_data$lof_z_mane, trans_data$lof_z_mane)
  y_max <- max(cis_data$lof_z_mane, trans_data$lof_z_mane)

  # Helper for violin border colors
  violin_colors <- cis_trans_color_scale
  violin_border_colors <- colorspace::darken(violin_colors, amount = 0.3)

  # Scatter plots
  remove_x_theme <- theme(
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

  scatter_cis <- ggplot(cis_data, aes(abs(beta), lof_z_mane, color = typeCisTrans)) +
    geom_point(shape = 20, size = 0.5, alpha = 0.2) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 0.3) +
    scale_colour_manual(values = cis_trans_color_scale, guide = NULL) +
    scale_x_log10(name = "|Beta|", limits = c(x_min, x_max)) +
    scale_y_continuous(name = "LoF Z-score", limits = c(y_min, y_max)) +
    geom_text(data = beta_lof_z_summary %>% filter(typeCisTrans == "cis"),
              aes(x = 1, y = 20, label = log_cor_label),
              parse = TRUE, size = 6 / .pt, inherit.aes = FALSE) +
    remove_x_theme +
    theme(panel.grid.major = element_line())

  scatter_trans <- ggplot(trans_data, aes(abs(beta), lof_z_mane, color = typeCisTrans)) +
    geom_point(shape = 20, size = 0.5, alpha = 0.2) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 0.3) +
    scale_colour_manual(values = cis_trans_color_scale, guide = NULL) +
    scale_x_log10(name = "|Beta|", limits = c(x_min, x_max)) +
    scale_y_continuous(name = "LoF Z-score", limits = c(y_min, y_max)) +
    geom_text(data = beta_lof_z_summary %>% filter(typeCisTrans == "trans"),
              aes(x = 1, y = 20, label = log_cor_label),
              parse = TRUE, size = 6 / .pt, inherit.aes = FALSE) +
    annotation_logticks(sides='b', linewidth=(0.5 / (ggplot2::.pt * 72.27/96)),
                        outside=T,
                        short = unit((0.1 / (ggplot2::.pt * 72.27/96)), "cm"),
                        mid = unit((0.1 / (ggplot2::.pt * 72.27/96)), "cm"),
                        long = unit((0.1 / (ggplot2::.pt * 72.27/96)), "cm"),) +
    coord_cartesian(clip="off") +
    remove_x_theme +
    theme(panel.grid.major = element_line())

  # Summary plot
  summary_plot <- ggplot(data %>% filter(any_CS),
                         aes(x = abs(beta),
                             y = ordered(typeCisTrans, levels = c("trans", "cis")),
                             fill = typeCisTrans,
                             color = typeCisTrans)) +
    geom_violin(width = 1, linewidth=0) +
    geom_boxplot(width = 0.05, notch = T, notchwidth = 0, outlier.shape = NA, fill = "black", color = "black") +
    scale_fill_manual(values = cis_trans_color_scale, guide = "none") +
    scale_color_manual(values = violin_border_colors, guide = "none") +
    scale_x_log10(name = "|Beta|", limits = c(x_min, x_max)) +
    ylab(NULL) + xlab(NULL) +
    theme(panel.grid.major = element_line())
    geom_text(data = beta_summary,
              aes(x = 1, y = 1, label = mwu_label),
              parse = TRUE, size = 6 / .pt, inherit.aes = FALSE)

  # Combine panels
  final_plot <- scatter_cis / scatter_trans / summary_plot +
    plot_layout(heights = c(2.5, 2.5, 1))

  return(final_plot)
}


rank_genes <- function(exp, tested_genes, method="RRA") {
  ranking_fun <- function(x){ranked <- rank(x, na.last = "keep")
    return(ranked)}

  # Include only genes that are detectable in 20% of cohorts
  exp <- exp[exp$gene %in% tested_genes,]

  ranks_genes <- apply(exp[, -1, with = FALSE], 2, ranking_fun)
  rownames(ranks_genes) <- exp$gene
  ranks_genes <- as.data.frame(ranks_genes)

  inp_list <- list()

  for (i in 1:ncol(ranks_genes)){
    temp <- data.table(Gene = rownames(ranks_genes), rank = ranks_genes[,i])
    temp <- temp[order(temp$rank, decreasing = TRUE),]
    temp <- temp[!is.na(temp$rank), ]

    temp_list <- list(temp$Gene)
    names(temp_list) <- colnames(ranks_genes)[i]

    inp_list <- c(inp_list, temp_list)

  }
  ranks <- aggregateRanks(inp_list, method=method, full = TRUE)

  all(order(ranks$Score) == 1:nrow(ranks))

  # Cut genes into 10 bins, based on expression
  genes_per_bin <- ceiling(nrow(ranks) / 10)
  last_bin <- nrow(ranks) - (genes_per_bin * 10)

  bins <- c(rep("Bin10", times = genes_per_bin),
            rep("Bin9", times = genes_per_bin),
            rep("Bin8", times = genes_per_bin),
            rep("Bin7", times = genes_per_bin),
            rep("Bin6", times = genes_per_bin),

            rep("Bin5", times = genes_per_bin),
            rep("Bin4", times = genes_per_bin),
            rep("Bin3", times = genes_per_bin),
            rep("Bin2", times = genes_per_bin),
            rep("Bin1", times = genes_per_bin + last_bin))

  ranks$bin <- ordered(bins, levels=paste0("Bin", 1:10))
  return(ranks)
}

# Declare function definitions
credible_set_locations_relative_to_gene_patchwork <- function(results_long_proc, window, gene_locus_colors) {

  # Setup
  upstream_breaks <- seq(window["start"], 0, window["by"])
  downstream_breaks <- seq(0, window["end"], window["by"])

  position_scales <- list(
    "upstream" = scale_x_continuous(breaks = upstream_breaks, labels = upstream_breaks / 1000),
    "intragenic" = scale_x_continuous(),
    "downstream" = scale_x_continuous(breaks = downstream_breaks, labels = downstream_breaks / 1000)
  )

  facet_names <- c(
    "upstream" = "Distance to TSS (kb)",
    "intragenic" = "Location within gene (%)",
    "downstream" = "Distance from TES (kb)"
  )

  regions <- c("upstream", "intragenic", "downstream")

  # Prepare summary data
  binned_data <- results_long_proc %>%
    filter(location_relative_to_gene != "trans") %>%
    group_by(location_relative_to_gene) %>%
    mutate(distance_capped_bins = cut(
      distance_capped, seq(0, max(abs(distance_capped)), by = max(abs(distance_capped))/20) * sign(max(distance_capped)),
      include.lowest = T)) %>%
    group_by(location_relative_to_gene, distance_capped_bins) %>%
      mutate(distance_capped_mean = mean(distance_capped))


    summary_data <- binned_data %>%
    group_by(location_relative_to_gene, distance_capped_bins) %>%
    summarise(
      n = n(),
      average_cs_rank = mean(as.numeric(as.character(cs_rank))),
      distance_capped = mean(distance_capped),
      .groups = "drop"
    )

  # Determine fixed y-axis limits
  max_histogram_count <- summary_data %>%
    pull(n) %>%
    max()

  max_summary_y <- summary_data %>%
    pull(average_cs_rank) %>%
    max() %>% ceiling()

  min_summary_y <- 1


  # Functions to create individual plots
  make_histogram_plot <- function(region, show_y_axis = FALSE) {
    y_axis_element <- if (show_y_axis) element_text() else element_blank()
    y_axis_line_element <- if (show_y_axis) element_line() else element_blank()
    y_axis_ticks_element <- if (show_y_axis) element_line() else element_blank()

    ggplot(
      results_long_proc %>% filter(location_relative_to_gene == region),
      aes(distance_capped, alpha = cs_rank, fill = location_relative_to_gene)
    ) +
      geom_histogram(position = "stack", bins=20, boundary=0) +
      scale_alpha_discrete(range = c(1, 0.1)) +
      scale_fill_manual(values = gene_locus_colors[region]) +
      position_scales[[region]] +
      labs(x = facet_names[region], y = if (show_y_axis) "Count" else NULL) +
      coord_cartesian(ylim = c(0, max_histogram_count)) +
      guides(fill = "none") +
      theme(
        legend.position = "right",
        panel.grid.major.y = element_line(color = "grey"),
        axis.line.y = y_axis_line_element,
        axis.title.y = y_axis_element,
        axis.text.y = y_axis_element,
        axis.ticks.y = y_axis_ticks_element
      )
  }

  make_summary_plot <- function(region, show_y_axis = FALSE) {
    y_axis_element <- if (show_y_axis) element_text() else element_blank()
    y_axis_ticks_element <- if (show_y_axis) element_line() else element_blank()
    y_axis_line_element <- if (show_y_axis) element_line() else element_blank()

    ggplot(
      binned_data %>% filter(location_relative_to_gene == region),
      aes(distance_capped_mean, as.numeric(as.character(cs_rank)), color = location_relative_to_gene)
    ) +
      stat_summary(fun.data = mean_se, fun.args = list(mult = 1), shape = 20,
                   size=(0.25 / (ggplot2::.pt * 72.27/96)), linewidth=(0.5 / (ggplot2::.pt * 72.27/96)),
                   geom = "pointrange") +
      geom_smooth(aes(color = -1), color = "#D55E00",
                  linewidth=(1 / (ggplot2::.pt * 72.27/96))) +
      scale_color_manual(values = gene_locus_colors[region]) +
      guides(color = "none") +
      position_scales[[region]] +
      labs(x = facet_names[region], y = if (show_y_axis) "Average CS rank" else NULL) +
      coord_cartesian(ylim = c(min_summary_y, max_summary_y)) +
      scale_y_continuous(breaks=c(1,2,3,4,5,6)) +
      theme(
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey"),
        axis.line.y = y_axis_line_element,
        axis.title.y = y_axis_element,
        axis.text.y = y_axis_element,
        axis.ticks.y = y_axis_ticks_element,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      )
  }

  # Generate all plots
  summary_plots <- mapply(make_summary_plot, regions, show_y_axis = c(TRUE, FALSE, FALSE), SIMPLIFY = FALSE)
  histogram_plots <- mapply(make_histogram_plot, regions, show_y_axis = c(TRUE, FALSE, FALSE), SIMPLIFY = FALSE)

  # Combine with patchwork
  combined_plot <- (wrap_plots(summary_plots, ncol = 3) / wrap_plots(histogram_plots, ncol = 3)) +
    plot_layout(guides = "collect", heights = c(1.2, 2)) &
    theme(legend.position = "right")

  return(combined_plot)
}

format_pvalues <- function(p.value, digits = 1, labelled = TRUE, alpha = 0.05) {
  sapply(p.value, function(p) {
    if (p > alpha) {
      return(ifelse(labelled, "~p==NS", "NS"))
    } else if (p == 0) {
      format_string <- paste0(ifelse(labelled, "~p<", "''<"), "%.", digits, "e")
      formatted <- sprintf(format_string, .Machine$double.xmin)
      gsub("e", " %*% 10^", formatted)
    } else {
      format_string <- paste0(ifelse(labelled, "~p==", ""), "%.", digits, "e")
      formatted <- sprintf(format_string, p)
      gsub("e", " %*% 10^", formatted)
    }
  })
}

variant_class_enrichments <- function(finemapped_variant_classes_unique_with_matched, classes, group) {

  variant_class_enrichment_mat <- finemapped_variant_classes_unique_with_matched %>%
    group_by(across(all_of(c(group, classes)))) %>%
    summarise(n = n()) %>%
    pivot_wider(id_cols = all_of(classes), names_from = group, values_from = "n") %>%
    rowwise() %>%
    filter(!is.na(yes) & !is.na(no)) %>%
    ungroup() %>%
    mutate(p_eqtl = yes / (no + yes) * 100,
           background_yes = sum(yes) - yes,
           background_no = sum(no) - no) %>%
    rowwise() %>%
    mutate(
      conf.mat = list(matrix(c(yes, no, background_yes, background_no), nrow = 2)),
      fischer = list(tidy(fisher.test(conf.mat)))
    ) %>%
    rename(setNames(c("variant_class"), c(classes))) %>%
    unnest(fischer, names_sep = "_") %>%
    mutate(variant_class_label = str_wrap(variant_class, width = 14)) %>%
    arrange(str_detect(variant_class, "hgsv")) %>%
    mutate(variant_class_label = ordered(variant_class_label, levels=variant_class_label))

  variant_class_enrichments <- ggplot(
    variant_class_enrichment_mat,
    aes(y = variant_class_label, x = fischer_estimate, xmin = fischer_conf.low, xmax = fischer_conf.high)) +
    geom_vline(xintercept = 1, color="grey70", size=(0.5 / (ggplot2::.pt * 72.27 / 96)), linetype="solid") +
    geom_col(orientation = "y", fill = "#C3E9FD") +
    geom_errorbar(width = 0.2, linewidth = (0.5 / (ggplot2::.pt * 72.27 / 96))) +
    geom_text(aes(label = format_pvalues(fischer_p.value, labelled = F), x = 0), parse = T, hjust = "inward", nudge_x = 0.08, size.unit = "pt", size = 6) +
    scale_y_discrete(name = "") +
    scale_x_continuous(name = "Enrichment (OR ±CI)")

  return(list("plot"=variant_class_enrichments, "data"=variant_class_enrichment_mat))
}


variant_class_enrichments_over_snps <- function(finemapped_variant_classes_unique_with_matched, col) {

  variant_class_enrichment_mat <- finemapped_variant_classes_unique_with_matched %>%
    group_by(across(all_of(c(col, "variant_class")))) %>%
    summarise(n = n()) %>%
    pivot_wider(id_cols = variant_class, names_from = col, values_from = "n") %>%
    rowwise() %>%
    filter(!is.na(yes) & !is.na(no)) %>%
    ungroup() %>%
    mutate(p_eqtl = yes / (no + yes) * 100,
           background = variant_class == "SNP",
           background_yes = yes[background],
           background_no = no[background]) %>%
    rowwise() %>%
    mutate(
      conf.mat = list(matrix(c(yes, no, background_yes, background_no), nrow = 2)),
      fischer = list(tidy(fisher.test(conf.mat)))
    ) %>%
    unnest(fischer, names_sep = "_") %>%
    mutate(variant_class_label = str_wrap(variant_class, width = 14)) %>%
    arrange(str_detect(variant_class, "hgsv")) %>%
    mutate(variant_class_label = ordered(variant_class_label, levels=variant_class_label))

    variant_class_enrichments <- ggplot(
    variant_class_enrichment_mat,
    aes(y = variant_class_label, x = fischer_estimate, xmin = fischer_conf.low, xmax = fischer_conf.high)) +
      geom_vline(xintercept = 1, color="grey70", size=(0.5 / (ggplot2::.pt * 72.27 / 96)), linetype="solid") +
      geom_col(orientation = "y", fill = "#C3E9FD") +
    geom_errorbar(width = 0.2, linewidth = (0.5 / (ggplot2::.pt * 72.27 / 96))) +
    geom_text(aes(label = format_pvalues(fischer_p.value, labelled = F), x = 0), parse = T, hjust = "inward", nudge_x = 0.08, size.unit = "pt", size = 6) +
    scale_y_discrete(name = "") +
    scale_x_continuous(name = "Enrichment (OR ±CI)")

  return(list("plot"=variant_class_enrichments, "data"=variant_class_enrichment_mat))
}

finemapping_per_variant <- function(per_variant_path, independent_variant_set_annot_filtered, variant_annotations) {
  value <- independent_variant_set_annot_filtered %>%
    filter(!is.na(pip)) %>%
    select(variant_index, phenotype, cluster, gene_cluster, CS_identifier, typeCisTrans)

  finemapping_per_variant_proc <- fread(per_variant_path) %>%
    select(variant_index, phenotype, cluster, gene_cluster, SusieRss_CS, SusieRss_pip, beta, standard_error) %>%
    left_join(value, by = join_by(variant_index, phenotype, cluster, gene_cluster)) %>%
    inner_join(variant_annotations, by = join_by(variant_index)) %>%
    group_by(chromosome, cluster, gene_cluster, SusieRss_CS, phenotype) %>%
    mutate(credible_set_retained = any(!is.na(CS_identifier)),
           lead_variant_in_cs = !is.na(CS_identifier),
           typeCisTrans = case_when(
             any(typeCisTrans == "cis") ~ "cis",
             any(typeCisTrans == "trans") ~ "trans",
             TRUE ~ ""
           )) %>%
    filter(credible_set_retained)

  return(finemapping_per_variant_proc)
}


# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Define variables
  outcome_variables <- c()

  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))
  location <- str_split(variables["available_genes"], "freeze3")[[1]][1]

  available_genes <- fread(variables["available_genes"], header=F)$V1
  outcome_variables["n_available_genes_meta_analysis"] <- length(available_genes)

  available_variants <- fread(variables["available_variants"], header=F)$V1

  # Window sizes
  cis_window <- 1e6
  twilight_window <- 5e6

  outcome_variables["cis_window"] <- cis_window
  outcome_variables["twilight_window"] <- twilight_window

  # Windows to plot
  window_size <- 1e6
  window <- c("start"=-window_size/2, "end"=window_size/2, "by" = 2e5)

  # Color pallete
  gene_locus_colors <- c("upstream" = "#D7817F", "intragenic" = "#D7817F", "downstream" = "#D7817F")
  cis_trans_color_scale <- c("cis" = "#CC79A7", "trans" = "#009E73")

  # Load variant annotations, has various annotations besides just positions, names and alleles
  variant_annotations <- read_parquet(variables["variant_annotations"])

  gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
    filter(type == "gene") %>% select(gene_id, gene_name, seqid, start, end, strand) %>%
    mutate(TSS = if_else(strand == "+", start, end))

  # Get calculated consequences
  calculated_consequences <- fread(
    variables["all_consequences"],
    header=F, col.names = c("term", "description", "accession", "display", "impact")) %>%
    mutate(display2 = str_remove(display, " variant"))

  ensembl_consequence_labels <- calculated_consequences %>% select(term, display2) %>% deframe()

  # Get background set of genes
  # I.e. the maximum set of genes for which we could hypothetically detect eQTLs
  # - The intersect of the genes in the gene-reference file
  # - and in the set of ensembl IDs for which the meta-analysis was run
  gene_reference_filtered <- gene_reference %>% filter(gene_id %in% available_genes)
  tested_genes <- gene_reference %>% filter(gene_id %in% available_genes) %>% pull(gene_id)

  # Filter the finemapped gene list on those that are not in the HLA or on CHR X
  finemapped_genes_cis <- gene_reference %>%
    filter(
      seqid != "X",
      gene_id %in% available_genes,
      !(seqid == 6 & start <= 34000000 & 25000000 <= end)
    ) %>%
    pull(gene_id)

  # Write a list of all tested genes
  fwrite(data.frame(gene = tested_genes), "tested_genes.txt", col.names=F, row.names=F, quote=F)

  # For trans, the regular set of tested genes can be used
  finemapped_genes_trans <- tested_genes

  # Count all tested genes
  outcome_variables["n_available_genes_meta_analysis_in_gene_reference"] <- sum(available_genes %in% gene_reference$gene_id)

  # Load the main file of independent variants
  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type))

  # Process all cis- and trans-eQTLs to put them in the supplement
  all_identified_eqtls_for_supplement <- independent_variant_set_annot_filtered %>%
    select(ensembl_gene_id=phenotype, gene_name, variant, chromosome, bp, eff_allele, non_eff_allele, beta, standard_error, p_value, pip, sample_size, type)

  # Write supplementary file
  fwrite(all_identified_eqtls_for_supplement, "all_identified_eQTLs_for_supplement_20251022.txt.gz", sep="\t", col.names=T, row.names=F)

  # Get all variants incluced in finemapped loci
  finemapping_per_variant_proc <- finemapping_per_variant(
    variables["finemapping_results"], independent_variant_set_annot_filtered, variant_annotations)

  # Write a set of credible set variants that pass filtering
  fwrite(finemapping_per_variant_proc %>% filter(typeCisTrans == "cis"), "cis-eQTLs_allCredibleSetVariants_OnlyKeptCredibleSets_20251028.txt.gz", sep="\t", quote=F, row.names=F, col.names=T)

  # Get a summary table that outlines for all CSs how often a credible set contains various variant classes,
  # and how often a variant class is the lead variant
  finemapping_variant_classes_per_cs <- finemapping_per_variant_proc %>%
    group_by(chromosome, cluster, gene_cluster, SusieRss_CS, phenotype, typeCisTrans) %>%
    summarise(any_indel = any(variant_class == "Insertion" | variant_class == "Deletion"),
              any_sv = any(str_detect(variant_class, "hgsv")),
              lead_variant_class = variant_class[lead_variant_in_cs],
              lead_variant_pip = SusieRss_pip[lead_variant_in_cs],
              lead_variant = variant[lead_variant_in_cs],
              all_pips = paste0(round(SusieRss_pip, digits=2), collapse=";"),
              all_classes = paste0(variant_class, collapse=";")
    )

  # Lead variant (finemapping-only) variant classes over all eQTLs
  finemapping_variant_classes_per_cs %>%
    group_by(lead_variant_class) %>%
    tally() %>%
    ungroup() %>%
    mutate(p = n / sum(n))

  # Lead variant (finemapping-only) variant classes over all unique variants
  finemapping_variant_classes_per_cs %>%
    ungroup() %>%
    distinct(lead_variant_class, lead_variant) %>%
    group_by(lead_variant_class) %>%
    tally() %>%
    ungroup() %>%
    mutate(p = n / sum(n))

  # Perform enrichment analyses to test if we see more indels/SNPs/SVs than what is expected based on all tested variants
  # Now load a set of matching cis and trans-eQTL variants
  matching_set_cis <- read_parquet(variables["matching_set_cis"]) %>%
    select(focal_index, variant_index) %>% mutate(matched_to_type="cis")
  matching_set_trans <- read_parquet(variables["matching_set_trans"]) %>%
    select(focal_index, variant_index) %>% mutate(matched_to_type="trans")

  # Concatenate all matching variants, but remove possible duplicate variants so that every eQTL variant is represented with an approximate equal number of control variants
  matching_set_any <- as.data.table(bind_rows(matching_set_cis, matching_set_trans %>% anti_join(matching_set_cis, by = c("focal_index")))) %>%
    inner_join(variant_annotations) %>% select(focal_index, variant, variant_class)

  # Get the set of unique focal variant ids
  matched_focal_indices <- unique(matching_set_any$focal_index)

  # For all focal variant ids, get the variant classes
  variant_class_of_focal_snps <- variant_annotations %>%
    filter(variant_index %in% matched_focal_indices) %>%
    select(variant, variant_index, variant_class) %>%
    mutate(variant_class_broad = if_else(str_detect(variant_class, "hgsv"), "SV", variant_class))

  # Merge the set of focal variant ids with variant classes to the eQTL data
  finemapped_variant_classes <- finemapping_variant_classes_per_cs %>%
    ungroup() %>%
    select(lead_variant, lead_variant_pip) %>%
    inner_join(
      variant_class_of_focal_snps, by = c("lead_variant" = "variant"))

  # First, test if there is an enrichment compared to control eQTLs. Also generates a figure
  finemapped_variant_classes_with_matched <- bind_rows(list(
    "yes" = finemapped_variant_classes,
    "no" = matching_set_any %>% filter(focal_index %in% finemapped_variant_classes$variant_index)), .id = "eqtl_variant")

  variant_class_enrichments_eqtls <- variant_class_enrichments(finemapped_variant_classes_with_matched, col="eqtl_variant")
  variant_class_enrichments_eqtls
  ggsave(sprintf("variant_class_enrichments_eqlts_%s.pdf", Sys.Date()), variant_class_enrichments_all$plot, width=50, height=50, units="mm")
  fwrite(variant_class_enrichments_all$data, sprintf("variant_class_enrichments_%s.txt", Sys.Date()), sep="\t", row.names=F, col.names=T)

  # Second, test if there is an enrichment compared to control variants. This is limited to just unique variants. Also generates a figure
  finemapped_variant_classes_unique <- finemapping_variant_classes_per_cs %>%
    ungroup() %>%
    select(lead_variant, lead_variant_pip) %>%
    distinct(lead_variant, .keep_all=T) %>%
    inner_join(
      variant_class_of_focal_snps, by = c("lead_variant" = "variant"))

  finemapped_variant_classes_unique_with_matched <- bind_rows(list(
    "yes" = finemapped_variant_classes_unique,
    "no" = matching_set_any %>% filter(focal_index %in% finemapped_variant_classes_unique$variant_index)), .id = "eqtl_variant") %>%
    mutate(variant_classes_very_broad = case_when(variant_class %in% c("Insertion", "Deletion") ~ "Indel", str_detect(variant_class, "(hgsv)") ~ "SV", variant_class == "SNP" ~ "SNP"))

  variant_class_enrichments_unique_variants <- variant_class_enrichments(finemapped_variant_classes_unique_with_matched, classes = "variant_classes_very_broad", group = "eqtl_variant")
  variant_class_enrichments_unique_variants
  ggsave(sprintf("variant_class_enrichments_unique_variants_%s.pdf", Sys.Date()), variant_class_enrichments_unique_variants$plot, width=50, height=50, units="mm")
  fwrite(variant_class_enrichments_unique_variants$data, sprintf("variant_class_enrichments_unique_variants_%s.txt", Sys.Date()), sep="\t", row.names=F, col.names=T)

  # Finally, test if there is an enrichment compared to all variants. This is limited to just unique variants. Also generates a figure
  background_unmatched <- variant_annotations %>% filter(variant_index %in% available_variants, !variant_index %in% independent_variant_set_annot_filtered$variant_index)

  finemapped_variant_classes_unique_with_unmatched_background <- bind_rows(list(
    "yes" = finemapped_variant_classes_unique,
    "no" = background_unmatched), .id = "eqtl_variant") %>%
    mutate(variant_classes_very_broad = case_when(variant_class %in% c("Insertion", "Deletion") ~ "Indel", str_detect(variant_class, "(hgsv)") ~ "SV", variant_class == "SNP" ~ "SNP"))

  variant_class_enrichments_all_unmatched <- variant_class_enrichments(finemapped_variant_classes_unique_with_unmatched_background, classes = "variant_classes_very_broad", group = "eqtl_variant")
  variant_class_enrichments_all_unmatched
  ggsave(sprintf("variant_class_enrichments_unmatched_background_%s.pdf", Sys.Date()), variant_class_enrichments_all_unmatched$plot, width=50, height=50, units="mm")
  fwrite(variant_class_enrichments_all_unmatched$data, sprintf("variant_class_enrichments_unmatched_background_%s.txt", Sys.Date()), sep="\t", row.names=F, col.names=T)

  # We can also do something similar confined to finemapped variants with a very high pip
  top_finemapped_variant_classes_unique <- finemapping_variant_classes_per_cs %>%
    ungroup() %>%
    filter(!is.na(lead_variant_pip) & lead_variant_pip > 0.9) %>%
    select(lead_variant, lead_variant_pip) %>%
    distinct(lead_variant, .keep_all=T) %>%
    inner_join(
      variant_class_of_focal_snps, by = c("lead_variant" = "variant"))

  top_finemapped_variant_classes_unique_with_matched <- bind_rows(list(
    "yes" = top_finemapped_variant_classes_unique,
    "no" = matching_set_any %>% filter(focal_index %in% top_finemapped_variant_classes_unique$variant_index)), .id = "eqtl_variant")

  variant_class_enrichments_top <- variant_class_enrichments(top_finemapped_variant_classes_unique_with_matched, col="eqtl_variant")
  variant_class_enrichments_top

  # Per eQTL type, list the number and percentages of different variant classes
  finemapping_variant_classes_per_cs_per_type <- finemapping_variant_classes_per_cs %>%
    group_by(typeCisTrans) %>%
    summarise(n=n(), n_any_indel = sum(any_indel), n_any_sv = sum(any_sv),
              lead_indel = sum((lead_variant_class == "Insertion" | lead_variant_class == "Deletion") & lead_variant_pip > 0.9),
              lead_sv = sum(str_detect(lead_variant_class, "hgsv") & lead_variant_pip > 0.9))

  print(finemapping_variant_classes_per_cs_per_type)

  finemapping_variant_classes_per_cs_per_type %>%
    ungroup() %>%
    summarise(across(where(is.numeric), sum),
              p_lead_sv = lead_sv / n * 100,
              p_lead_indel = lead_indel / n * 100)


  # Estimate how much larger beta's are of cis-eQTLs compared to trans-eQTls
  # First, get just the lead cis- and trans-effects per gene
  lead_effects <- independent_variant_set_annot_filtered %>%
    group_by(phenotype, typeCisTrans) %>%
    slice_max(abs(beta/standard_error), with_ties=F) %>%
    ungroup()

  # Get metrics of averages for the absolute beta's
  lead_effects %>% group_by(typeCisTrans) %>%
    summarise(tidy(summary(abs(beta))))

  # Estimate if there is a difference
  res <- wilcox.test(abs(lead_effects$beta[lead_effects$typeCisTrans=="cis"]),
                     abs(lead_effects$beta[lead_effects$typeCisTrans=="trans"]))

  # Reformat the set of eQTLs:
  # Primarily, add the location of the gene relative to the variant
  # Order CSs by significance
  results_long_proc <- independent_variant_set_annot_filtered %>%
    mutate(
      typeCisTrans = if_else(type == "intermediate", "trans", type),
      z_score = beta/standard_error,
      bp = if_else(strand == "+", bp, -bp),
      start_ = if_else(strand == "+", start, -end),
      end = if_else(strand == "+", end, -start),
      start = start_) %>%
    mutate(
      distance_to_tss = if_else(chromosome == seqid, bp - start, NA_integer_),
      distance_to_tes = if_else(chromosome == seqid, bp - end, NA_integer_),
      location_relative_to_gene = ordered(case_when(
        chromosome != seqid ~ "trans",
        distance_to_tss < -5e6 ~ "trans",
        distance_to_tss > 5e6 ~ "trans",
        distance_to_tss < 0 ~ "upstream",
        distance_to_tes > 0 ~ "downstream",
        TRUE ~ "intragenic"
      ), levels=c("upstream", "intragenic", "downstream", "trans")),
      location_relative_to_gene_fine = ordered(case_when(
        distance_to_tss <= -1e7 ~ "trans",
        distance_to_tss <= -1e6 ~ "1-5mb upstream",
        distance_to_tss <= -5e4 ~ "50kb-1mb upstream",
        distance_to_tss <= 0 ~ "<50 kb upstream",
        distance_to_tes >= 1e7 ~ "trans",
        distance_to_tes >= 1e6 ~ "1-5mb downstream",
        distance_to_tes >= 5e4 ~ "50kb-1mb downstream",
        distance_to_tes >= 0 ~ "<50 kb downstream",
        chromosome != seqid ~ "trans",
        TRUE ~ "intragenic"
      ), levels=c("1-5mb upstream", "50kb-1mb upstream", "<50 kb upstream", "intragenic", "<50 kb downstream", "50kb-1mb downstream", "1-5mb downstream", "trans")),
      distance = case_when(
        location_relative_to_gene == "trans" ~ NA_integer_,
        location_relative_to_gene == "upstream" ~ distance_to_tss,
        location_relative_to_gene == "downstream" ~ distance_to_tes,
        location_relative_to_gene == "intragenic" ~ (bp - start) / (end - start) * 100
      ),
      distance_capped = pmax(pmin(distance, window["end"]), window["start"])) %>%
    group_by(phenotype, seqid) %>%
    mutate(cs_rank = pmin(row_number(dplyr::desc(abs(beta/standard_error))), 10)) %>%
    ungroup() %>%
    mutate(cs_rank = ordered(cs_rank))

  # For each coarse location of the variant relative to the gene, summarise the distribution fo CS ranks
  results_long_proc %>%
    mutate(cs_rank = as.numeric(cs_rank)) %>%
    group_by(location_relative_to_gene_fine) %>%
    summarise(min(cs_rank), max(cs_rank), mean(cs_rank), median(cs_rank), n())

  # CS rank higher upstream, intragenic, or downstream
  pairwise_mwu_res <- stats::pairwise.wilcox.test(
      as.numeric(results_long_proc$cs_rank), results_long_proc$location_relative_to_gene_fine,
      p.adjust.method = "bonferroni"
  )

  # Relationship between credible set rank and distance to gene body
  results_long_proc %>% filter(typeCisTrans == "cis", location_relative_to_gene != "intragenic") %>%
    select(cs_rank, distance, location_relative_to_gene) %>%
    group_by(location_relative_to_gene) %>%
    nest() %>%
    summarise(rank_dist_cor = map(
      data, ~ tidy(cor.test(as.numeric(as.character(.x$cs_rank)), abs(.x$distance)), method="spearman"))) %>%
    unnest(rank_dist_cor) %>%
    mutate(
      cor_label =
        paste0(sprintf("rho==%.2f", estimate), "*", format_pvalues(p.value)))

  # Get average distance to TSS and TSS
  results_long_proc %>%
    filter(typeCisTrans == "cis") %>%
    summarise(median(abs(distance_to_tss)),
              median(abs(distance_to_tes)))

  # Calculate how many cis eQTLs are in each location bin
  results_long_proc %>%
    filter(typeCisTrans == "cis") %>%
    group_by(location_relative_to_gene_fine) %>%
    tally() %>% ungroup() %>% mutate(percent = n / sum(n) * 100)

  # Plot the location of the variants relative to the linked genes for cis-eQTLs
  fineampping_cs_ranks_plot <- credible_set_locations_relative_to_gene_patchwork(results_long_proc, window, gene_locus_colors)
  # Combined plot
  ggsave(sprintf("finemapping_cs_ranks_%s.pdf", Sys.Date()), fineampping_cs_ranks_plot, width=104, height=58, units="mm")

  gene_ranks <- fread(variables["expression_bins"])

  # Get gnomad results
  # Use the mane gene transcript if available,
  # otherwise, use the score that indicates the largest gene conservation
  gnomad4 <- read.delim(variables["gnomad_constraint_metrics"])
  geneMetrics4 <- gnomad4 %>%
    group_by(gene_id) %>%
    mutate(lof_z = max(lof.z_score, na.rm=T),
           pLI = max(lof.pLI, na.rm=T),
           mis_z = max(mis.z_score, na.rm=T),
           syn_z = max(syn.z_score, na.rm=T),
           loeuf = min(lof.oe_ci.upper, na.rm=T),
           mane_select = mane_select == "true"
    ) %>%
    slice_max(mane_select, n = 1, with_ties = F) %>%
    mutate(loeuf_mane = if_else(mane_select, lof.oe_ci.upper, loeuf),
           pLI_mane = if_else(mane_select, lof.pLI, pLI),
           lof_z_mane = if_else(mane_select & !is.na(lof.z_score), lof.z_score, lof_z)) %>%
    select(gene_id, lof_z, lof_z_mane, mis_z, syn_z, loeuf, loeuf_mane, pLI, pLI_mane, mane_select) %>%
    mutate(across(c("lof_z", "pLI", "mis_z", "syn_z"), ~ if_else(.x == -Inf, NA, .x) )) %>%
    mutate(across(c("loeuf"), ~ if_else(.x == Inf, NA, .x) ))

  fwrite(geneMetrics4, "/Users/cawarmerdam/Documents/projects/eQTLGen/processed_data/gnomad_constraint_metrics_summarised_per_gene_id.tsv", sep="\t", row.names=F, col.names=T)

  # Merge gene metrics with expression ranks
  gene_metrics_proc_4 <- gene_ranks %>%
    left_join(geneMetrics4, by = c("Name" = "gene_id")) %>%
    pivot_longer(c(Score_cis, Score_trans, bin_cis, bin_trans), names_sep="_", names_to=c(".value", "typeCisTrans")) %>%
    filter(!is.na(bin))

  # Per finemapping locus, count the number of credible sets
  per_gene_finemapping_locus_results <- results_long_proc %>%
    mutate(
      typeCisTrans = if_else(type == "intermediate", "trans", type),
    ) %>%
    group_by(locus_str) %>%
    summarise(n_CS_SuSiE = sum(!is.na(pip)),
              n_CS = n()) %>%
    ungroup() %>%
    summarise(
      n=n(),
      n_any_CS_SuSiE = sum(n_CS_SuSiE > 0),
      n_over_1_CS_SuSiE = sum(n_CS_SuSiE > 1),
      n_CS_total_SuSiE = sum(n_CS_SuSiE),
      n_CS_total = sum(n_CS),
      mean_n_CS = mean(n_CS_SuSiE))

  print(per_gene_finemapping_locus_results)

  # Per gene, eQTL type, count the number of CSs (confined to SuSiE or not)
  per_gene_type_combination <- results_long_proc %>%
    mutate(
      varianceExplained = 2 * eaf * (1 - eaf) * (beta * beta),
      typeCisTrans = if_else(type == "intermediate", "trans", type),
      beta_bin = ordered(
        ntile(abs(beta), 5),
        labels=c("0-20th perc.", "20-40th perc.", "40-60th perc.", "60-80th perc.", "80-100th perc.")),
    ) %>%
    group_by(phenotype, typeCisTrans) %>%
    mutate(n_CSs_per_gene = n(),
           n_CSs_per_gene_SuSiE = sum(!is.na(pip)),
           varianceExplainedSum = min(sum(varianceExplained), 1),
           n_loci = n_distinct(finemapping_index)) %>%
    slice_max(abs(beta/standard_error), with_ties = F) %>%
    ungroup() %>%
    right_join(gene_metrics_proc_4, by=c("phenotype" = "Name", "typeCisTrans" = "typeCisTrans")) %>%
    mutate(
      n_CSs_per_gene = if_else(is.na(n_CSs_per_gene), 0, n_CSs_per_gene),
      n_CSs_per_gene_capped = if_else(n_CSs_per_gene > 10, 11, n_CSs_per_gene),
      n_CSs_per_gene_ord_capped = ordered(n_CSs_per_gene_capped, labels=c(0:10, ">10")),
      any_CS = n_CSs_per_gene > 0,
      n_loci = if_else(is.na(n_loci), NA, n_loci),
      n_CSs_per_gene_SuSiE = if_else(is.na(n_CSs_per_gene_SuSiE), NA, n_CSs_per_gene_SuSiE),
      any_CS_SuSiE = n_CSs_per_gene_SuSiE > 0) %>%
    group_by(bin, typeCisTrans) %>%
      mutate(ve_bin = ordered(
        ntile(abs(varianceExplainedSum), 5),
        labels=c("0-20th perc.", "20-40th perc.", "40-60th perc.", "60-80th perc.", "80-100th perc.")))

  # Now, calculate an average across the cis and trans eQTL types
  summary_cis_trans_finemapping <- per_gene_type_combination %>%
    group_by(typeCisTrans) %>%
    summarise(n_genes = n(),
              n_finemapping_loci = sum(n_loci, na.rm=T),
              mean_finemapping_loci = mean(n_loci, na.rm=T),
              mean_n_SCS = mean(n_CSs_per_gene_SuSiE, na.rm=T),
              mean_n_CS = mean(n_CSs_per_gene, na.rm=T),
              n_over_1_CS = sum(n_CSs_per_gene_SuSiE > 1, na.rm=T),
              p_over_1_CS_genes = n_over_1_CS/n_genes * 100,
              p_over_1_CS_loci = n_over_1_CS/n_finemapping_loci * 100,
              n_any_CS = sum(any_CS),
              p_any_CS = n_any_CS / n_genes * 100,
              n_any_CS_SuSiE = sum(any_CS_SuSiE, na.rm=T),
              p_any_CS_SuSiE = n_any_CS_SuSiE / n_genes * 100,
              n_total_CS_SuSiE = sum(n_CSs_per_gene_SuSiE, na.rm=T),
              n_total_CS = sum(n_CSs_per_gene))

  print(summary_cis_trans_finemapping, width=Inf)

  # Also list CS sizes
  results_CS_sizes <- results_long_proc %>%
    group_by(is.na(pip), typeCisTrans) %>%
    summarise(n_CSs = n(),
              n_CS_size_any = sum(!is.na(CS_size)),
              n_CSS = sum(!is.na(pip)),
              mean_CS_size = mean(CS_size),
              n_CS_size_exactly_1 = sum(CS_size == 1),
              n_CS_size_exactly_1_pip_over_90p = sum(CS_size == 1 & pip > 0.9))

  print(results_CS_sizes, width=Inf)

  # Write gene summary for the supplement
  per_gene_type_combination_for_supplement <- per_gene_type_combination %>%
    group_by(phenotype) %>%
    mutate(gene_name = case_when(any(!is.na(gene_name)) ~ na.omit(gene_name)[1], TRUE ~ gene_name)) %>%
    ungroup() %>%
    select(phenotype, gene_name, typeCisTrans, n_CSs_per_gene, bin, Score) %>%
    pivot_wider(values_from=c(n_CSs_per_gene, bin, Score), names_from="typeCisTrans") %>%
    left_join(geneMetrics4, by = c("phenotype" = "gene_id"))

  fwrite(per_gene_type_combination_for_supplement, "per_gene_type_combination_for_supplement.txt", sep="\t", col.names=T, row.names=F)

  # Plot the relationship between the number of eQTL credible sets, expression and gene conservation
  eqtls_per_gene <- eqtl_overview_plot(per_gene_type_combination)
  eqtls_per_gene
  ggsave("eqtls_per_gene_v2.pdf", eqtls_per_gene, width=130, height=86, units="mm")

  # Plot the relationship between the number of eQTL credible sets, expression and gene conservation,
  # but using violin plots, and for supplements
  p_sup_cis <- eqtl_credible_set_leouf_supplement(per_gene_type_combination %>% filter(typeCisTrans == "cis"))
  p_sup_trans <- eqtl_credible_set_leouf_supplement(per_gene_type_combination %>% filter(typeCisTrans == "trans"))

  ggsave("eqtls_per_gene_supplement_cis.png", p_sup_cis, width=180, height=180, units="mm")
  ggsave("eqtls_per_gene_supplement_trans.png", p_sup_trans, width=180, height=180, units="mm")

  # Perform tests to understand the relationships were are seeing.
  per_gene_type_combination_cis <- per_gene_type_combination %>% filter(typeCisTrans == "cis") %>%
    filter(!is.na(lof_z_mane), is.finite(lof_z_mane)) %>%
    select(phenotype, any_CS, lof_z_mane, Score, bin, pLI_mane) %>%
    mutate(bin_num = as.numeric(bin))

  # Check if there is a difference in expression bin for the fine-mapped eQTLs and the eQTLs that were not fine-mapped
  wilcox.test(per_gene_type_combination_cis$bin_num[per_gene_type_combination_cis$any_CS], per_gene_type_combination_cis$bin_num[!per_gene_type_combination_cis$any_CS])
  ggplot(per_gene_type_combination_cis, aes(any_CS, Score)) + geom_violin() + geom_boxplot(width=0.1)

  # Check if there is a difference for the gene expression bins with the bins 6-10
  per_gene_type_combination_cis_bin6_10 <- per_gene_type_combination_cis %>% filter(bin %in% paste0("Bin",6:10))

  wilcox.test(per_gene_type_combination_cis$pLI_mane[per_gene_type_combination_cis$any_CS], per_gene_type_combination_cis$pLI_mane[!per_gene_type_combination_cis$any_CS])
  wilcox.test(per_gene_type_combination_cis_bin6_10$pLI_mane[per_gene_type_combination_cis_bin6_10$any_CS], per_gene_type_combination_cis_bin6_10$pLI_mane[!per_gene_type_combination_cis_bin6_10$any_CS])
  ggplot(per_gene_type_combination_cis, aes(any_CS, pLI_mane)) + geom_violin() + geom_boxplot(width=0.1)

  # Check if the relationship is different if we look at variants that are not colocalizing with blood cell type.
  downloads <- fread(variables["variant_coloc_classes"]) %>%
    filter(type %in% c("trans", "intermediate")) %>%
    group_by(phenotype) %>%
    summarise(n_ctc = sum(variant_type  == "Some blood-cell"), all_eqtls = n())

  # Check the relationship for the most expressed genes. This is the rightmost panel in figure 2b
  per_gene_type_combination_proc_bin10_cis <- per_gene_type_combination %>% filter(typeCisTrans == "cis", bin == "Bin10") %>%
    left_join(downloads, by = c("phenotype")) %>%
    filter(!is.na(loeuf_mane), is.finite(loeuf_mane))

  # Plot the rightmost panel in figure 2b
  bin10_loeuf <- ggplot(per_gene_type_combination_proc_bin10_cis,
         aes(n_CSs_per_gene, loeuf_mane, fill = n_CSs_per_gene_ord_capped, group = n_CSs_per_gene_ord_capped)) +
    geom_violin(linewidth=0) +
    geom_boxplot(width=0.1, fill="white") +
    scale_fill_manual(name = "Number of eQTLs", values = c("grey90", color_func(11)), guide=NULL) +
    scale_x_continuous(name = "number of eQTLs") +
    coord_cartesian(clip = "on") +
    scale_y_continuous(name = "LOEUF", breaks = c(0, 0.5, 1.0, 1.5, 2)) +
    geom_smooth(aes(n_CSs_per_gene, loeuf_mane), method = "lm", color = "#D55E00", inherit.aes=F,
                linewidth=(1 / (ggplot2::.pt * 72.27/96)))

  ggsave("n_eqtls_vs_loeuf_bin10_cis.pdf", bin10_loeuf, width=72, height=54, units="mm")

  # Check what is happening for trans bin4
  per_gene_type_combination_proc_bin4_trans <- per_gene_type_combination %>% filter(typeCisTrans == "trans", bin == "Bin4") %>%
    left_join(downloads, by = c("phenotype")) %>%
    mutate(n_ctc = replace_na(n_ctc, 0), all_eqtls = replace_na(all_eqtls, 0),
           n_non_ctc = all_eqtls - n_ctc) %>%
    filter(!is.na(loeuf_mane), is.finite(loeuf_mane))

  # Plot the relationshop for trans bin4
    ggplot(per_gene_type_combination_proc_bin4_trans,
           aes(n_CSs_per_gene, loeuf_mane, color = n_CSs_per_gene_ord_capped)) +
    geom_point(color = "black", alpha = 0.1) +
    stat_summary(fun.data = mean_se, fun.args = list(mult = 1), shape = 20,
    size=(0.25 / (ggplot2::.pt * 72.27/96)), linewidth=(0.5 / (ggplot2::.pt * 72.27/96)),
    geom = "pointrange") +
    scale_x_continuous(name = "number of eQTLs") +
    coord_cartesian(clip = "on") +
    scale_y_continuous(name = "LOEUF", breaks = c(0.5, 0.75, 1.0, 1.25)) +
    geom_smooth(method = "lm", aes(color = -1), color = "#D55E00",
    linewidth=(1 / (ggplot2::.pt * 72.27/96)))

  per_gene_type_combination_bin4_trans <- per_gene_type_combination %>% filter(typeCisTrans == "trans", bin == "Bin4")

  cor.test(per_gene_type_combination_proc_bin4_trans$n_CSs_per_gene, per_gene_type_combination_proc_bin4_trans$loeuf_mane, method = "spearman")
  per_gene_type_combination_bin4_trans_within_10_eqtls <- per_gene_type_combination_bin4_trans %>%
    filter(n_CSs_per_gene < 11)

  cor.test(per_gene_type_combination_bin4_trans_within_10_eqtls$n_CSs_per_gene, per_gene_type_combination_bin4_trans_within_10_eqtls$loeuf_mane, method = "spearman")

  res <- per_gene_type_combination_bin4_trans_within_10_eqtls %>%
    group_by(n_CSs_per_gene) %>%
    summarise(loeuf_mane_mean = mean(loeuf_mane, na.rm=T))

  cor.test(res$n_CSs_per_gene, res$loeuf_mane_mean, method = "spearman")

  # Plot the relationship between beta's, LoF Z-scores and eQTL type
  beta_lof_z_plot <- plot_beta_lof_summary_split_facets(
    data = per_gene_type_combination %>% ungroup(),
    cis_trans_color_scale = cis_trans_color_scale
  )

  print(beta_lof_z_plot)

  ggsave("beta_to_lof_z.pdf", beta_lof_z_plot, width=48, height=72, units="mm")

  per_gene_type_combination %>% pivot_wider(id_cols = c(phenotype, bin), names_from=typeCisTrans, values_from = n_CSs_per_gene) %>%
  ggplot(aes(trans, cis, color=bin)) +
    geom_jitter() +
    facet_wrap(vars(bin))

  per_gene_type_combination %>%
    group_by(typeCisTrans, n_CSs_per_gene) %>%
    summarise(n_genes = n()) %>%
  ggplot(aes(n_CSs_per_gene, n_genes, color=typeCisTrans)) +
    geom_line(aes(group=typeCisTrans), color="grey") +
    geom_point(shape=20) +
    scale_colour_manual(values= cis_trans_color_scale) +
    scale_x_continuous(name="Number of eQTLs per gene", trans = "pseudo_log", breaks=c(0,1,5,10,20,40,60)) +
    scale_y_continuous(name="N genes") +
    guides(x="axis_logticks") +
    coord_cartesian(clip="off") +
    theme(legend.position = "none")

  ggsave("top_eqtls_per_gene.pdf", top_eqtl_genes, width=50, height=80, units="mm")

  top_cis_trans <- per_gene_type_combination %>% group_by(typeCisTrans) %>%
    mutate(top_cis_trans_rank = row_number(dplyr::desc(n_CSs_per_gene))) %>%
    ungroup() %>%
    group_by(phenotype, gene_name)

  top_trans <- top_cis_trans %>% filter(any(top_cis_trans_rank < 20 & typeCisTrans=="trans")) %>%
    ungroup() %>%
    arrange(top_cis_trans_rank) %>%
    mutate(gene_label = if_else(!is.na(gene_name), gene_name, phenotype),
           gene_label = ordered(gene_label, levels=rev(unique(gene_label)))) %>%
    select(typeCisTrans, phenotype, gene_name, gene_label, n_CSs_per_gene)

  top_eqtl_genes <- ggplot(top_trans, aes(gene_label, n_CSs_per_gene, fill=typeCisTrans)) +
    geom_col() +
    scale_y_continuous(name="Number of eQTLs per gene") +
    scale_x_discrete(name="") +
    scale_fill_manual(name="",values=cis_trans_color_scale) +
    coord_flip() +
    theme(legend.position = "bottom")

  top_eqtl_genes

  ggsave("top_eqtls_per_gene.pdf", top_eqtl_genes, width=50, height=80, units="mm")

  per_gene_type_combination_summarised <- per_gene_type_combination %>%
    group_by(bin, n_CSs_per_gene_ord_capped, typeCisTrans) %>%
    summarise(n = n(),
              n_non_na_z = sum(!is.na(lof.z_score)),
              n_non_na_pli = sum(!is.na(lof.pLI)),
              median_pli = median(lof.pLI, na.rm=T),
              avg_lof.z_score = mean(lof.z_score, na.rm=T),
              se_lof.z_score = sd(lof.z_score, na.rm=T)/sqrt(n_non_na_z),
              avg_lof.z_score_lower = avg_lof.z_score - se_lof.z_score,
              avg_lof.z_score_upper = avg_lof.z_score + se_lof.z_score)

  # Define table with all cis-eQTLs and their annotations to match on
  # Define table with all background cis-eQTLs and their annotations to match on

  # Define table with all trans-eQTls and their annotations to match on
  # Define table with all background cis-eQTLs and their annotations to match on

    annotated_finemapped_variants <- left_join(
      variant_annotations,
      independent_variant_set_annot_filtered %>%
        distinct(variant_index, typeCisTrans))

  # Here calculate VIP enrichments
  matching_focal_variants <- unique(c(matching_set_cis$focal_index, matching_set_trans$focal_index))

  vip_largest_cons_raw <- fread(variables["vip_consequences"], select=c("ID", "Consequence")) %>%
    select(variant = ID, consequence = Consequence)

  vip_largest_cons <- vip_largest_cons_raw %>%
    separate_rows(consequence, sep = "&") %>%
    mutate(consequence_as_ordered = ordered(consequence, levels=rev(names(ensembl_consequence_labels)))) %>%
    group_by(variant) %>%
    slice_max(consequence_as_ordered, with_ties = F, n=1) %>%
    right_join(variant_annotations, by = c("variant")) %>%
      mutate(consequence = ifelse(is.na(consequence), "intergenic_variant", consequence))

  variant_consequnce_per_eqtl <- independent_variant_set_annot_filtered %>%
    inner_join(vip_largest_cons_raw) %>%
    mutate(largest_predicted_consequence = ifelse(is.na(consequence), "intergenic_variant", consequence)) %>%
    select(phenotype, variant, variant_index, type, largest_predicted_consequence)

  fwrite(variant_consequnce_per_eqtl, "largest_predicted_consequence_finemapped_eqtls_20251020.tsv.gz", sep="\t", row.names=F, col.names=T)

  vip_actual <- annotated_finemapped_variants %>%
    filter(!is.na(typeCisTrans)) %>%
    distinct(variant_index, typeCisTrans) %>%
    inner_join(vip_largest_cons) %>%
    mutate(matched = variant_index %in% matching_focal_variants)

  vip_matched_trans <- matching_set_trans %>% select(variant_index) %>%
    inner_join(vip_largest_cons) %>%
    mutate(typeCisTrans = "trans")

  vip_matched_cis <- matching_set_cis %>% select(variant_index) %>%
    inner_join(vip_largest_cons) %>%
    mutate(typeCisTrans = "cis")

  vip_matched <- bind_rows(vip_matched_cis, vip_matched_trans)

  vip_matched_control_dataset <- bind_rows(list("actual"=vip_actual, "matched"=vip_matched), .id="group") %>%
    group_by(consequence, typeCisTrans, group) %>%
    summarise(
      n_observed = n(),
      n_with_matched = sum(is.na(matched) | matched))

  fwrite(vip_matched_control_dataset, "variant_consequence_counts_per_group_20251020.txt.gz", sep="\t", col.names=T, row.names=F)
  vip_matched_control_dataset <- fread("~/Documents/projects/eQTLGen/users/umcg-rwarmerdam/variant_consequence_counts_per_group_20250612.txt.gz")

  fisher_enrichments <- vip_matched_control_dataset %>%
    pivot_wider(id_cols = c(consequence, typeCisTrans), names_from = group, values_from = c(n_with_matched, n_observed), values_fill = 0) %>%
    group_by(typeCisTrans) %>%
    mutate(
      background_actual = sum(n_with_matched_actual) - n_with_matched_actual,
      background_matched = sum(n_with_matched_matched) - n_with_matched_matched) %>%
    rowwise() %>%
    mutate(
      conf.mat = list(matrix(c(n_with_matched_actual, n_with_matched_matched, background_actual, background_matched), nrow = 2)),
      fischer = list(tidy(fisher.test(conf.mat)))
    ) %>%
    unnest(fischer, names_sep = "_") %>%
    ungroup() %>%
    arrange(dplyr::desc(fischer_estimate)) %>%
    mutate(
      consequence = recode(consequence, !!!ensembl_consequence_labels),
      consequence = ordered(
        consequence,
        levels = unique(c(consequence, "Intergenic"), fromLast = F)))

  fwrite(fisher_enrichments, "variant_consequence_enrichments_20251020.txt.gz", sep="\t", col.names=T, row.names=F)


  labelled_consequences <- fisher_enrichments %>%
    select(consequence, typeCisTrans, n_observed_actual) %>%
    filter(typeCisTrans %in% c("cis", "trans")) %>%
    pivot_wider(names_from = typeCisTrans, values_from = n_observed_actual, values_fill = 0) %>%
    mutate(
      consequence_labeled = paste0(consequence, "\n(cis: ", format(cis, big.mark = ",", justify="none", trim=T), "; trans: ", format(trans, big.mark = ",", trim=T), ")")
    )

  # Join back to main table and replace consequence
  fisher_enrichments2 <- fisher_enrichments %>%
    left_join(labelled_consequences, by = "consequence") %>%
    mutate(consequence2 = ordered(consequence, levels = unique(consequence)),
           n_tests = sum(n_with_matched_actual > 0 & consequence != "Intergenic"),
           bonferroni_alpha = 0.05 / n_tests,
           formatted_p = if_else(fischer_p.value < bonferroni_alpha, "'*'", ""))

  fwrite(fisher_enrichments2, "variant_consequence_enrichments_20251020.txt.gz", sep="\t", col.names=T, row.names=F)


  fisher_enrichments %>% group_by(typeCisTrans) %>%
    summarise(n())

  # Create the table-like plot
  table_plot <- ggplot(fisher_enrichments2 %>% group_by(consequence2) %>% filter(any(n_with_matched_actual > 10)),
                       aes(x=1, y = consequence2, colour=typeCisTrans)) +
    geom_text(aes(label = format(n_observed_actual, big.mark = ",")),
              hjust = 0, position=position_dodge2(preserve = "single", width=1), size=4, size.unit = "pt", parse=F) +
    scale_x_continuous(limits = c(0.9, 1.1), breaks = c(), name="", expand=c(0,0)) +
    scale_y_discrete(name = "") +
    scale_colour_manual(values=cis_trans_color_scale) +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      axis.text.x = element_blank(),
      legend.position = "none"
    )

  # Create the table-like plot
  table_plot_sign <- ggplot(fisher_enrichments2 %>% group_by(consequence2) %>% filter(any(n_with_matched_actual > 10)),
                       aes(x=1, y = consequence2, colour=typeCisTrans)) +
    geom_text(aes(label = formatted_p),
              hjust = 0, position=position_dodge2(preserve = "single", width=1), size=4, size.unit = "pt", parse=T) +
    scale_x_continuous(limits = c(0.9, 1.1), breaks = c(), name="", expand=c(0,0)) +
    scale_y_discrete(name = "") +
    scale_colour_manual(values=cis_trans_color_scale) +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      axis.text = element_blank(),
      legend.position = "none",
      axis = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank()
    )

  # Enrichments
  variant_consequence_enrichments <- ggplot(
    fisher_enrichments2 %>% filter(n_with_matched_actual > 10),
    aes(y = consequence2, x = fischer_estimate, xmin = fischer_conf.low, xmax = fischer_conf.high, fill=typeCisTrans)) +
    geom_vline(xintercept = 1, color="grey70", linewidth=(0.5 / (ggplot2::.pt * 72.27 / 96)), linetype="solid") +
    geom_col(orientation = "y", position=position_dodge2(preserve = "single")) +
    scale_fill_manual(values=cis_trans_color_scale) +
    geom_errorbar(width = 0.2, linewidth = (0.5 / (ggplot2::.pt * 72.27 / 96)), position=position_dodge(width=.9)) +
    scale_y_discrete(name = "") +
    scale_x_continuous(name = "Enrichment over intergenic\n(OR ±CI)", transform = pseudo_log_trans(base=10), limits=c(0, NA), breaks=c(0, 1, 2, 5, 10, 15)) +
    coord_cartesian(xlim=c(1,16)) +
    theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "null"), plot.background = element_blank())

  combined_plot <- table_plot + table_plot_sign + variant_consequence_enrichments +
    plot_layout(widths = c(1, 0.2, 6))

  print(combined_plot)
  ggsave("variant_consequences_with_raw_n.pdf", combined_plot, width=80, height=112, units="mm")

  variant_consequence_enrichments

  ggsave("variant_consequences.pdf", width=72, height=112, units="mm")

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}