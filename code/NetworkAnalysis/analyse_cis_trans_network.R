#!/usr/bin/env Rscript
#' Analysis of cis-trans network:
#' Specifically: Drawing a distribution plot of in-degree and out-degree
#' And checking enrichments of GWAS colocalization among trans-eQTLs that do have cis-effects vs those that do not.

# Load libraries
library(tidyverse)
library(data.table)
library(igraph)
library(poweRlaw)
library(broom)


# Declare constants
old <- theme_set(theme_classic(base_size = 6))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  text = element_text(family="Helvetica"),
  title = element_text(colour = "#595A5C", face = "plain")
)


cc_traits <- c("UKB_460K.blood_PLATELET_COUNT", "Eosinophil_count", "Eosinophil_percentage_of_white_cells", "Lymphocyte_percentage_of_white_cells",
               "Neutrophil_percentage_of_white_cells", "UKB_460K.blood_EOSINOPHIL_COUNT", "Neutrophil_count", "UKB_460K.blood_WHITE_COUNT",
               "White_blood_cell_count", "Monocyte_percentage_of_white_cells", "Hematocrit", "Hemoglobin", "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
               "UKB_460K.blood_RED_COUNT", "Basophil_count", "Basphil_percentage_of_leucocytes", "High_light_scatter_reticulocyte_percentage_of_red_cells",
               "Immature_fraction_of_reticulocytes", "Mean_corpuscular_hemoglobin", "Mean_corpuscular_volume", "Mean_reticulocyte_volume", "Mean_spheric_corpuscular_volume",
               "Plateletcrit", "Mean_corpuscular_hemoglobin_concentration", "Mean_platelet_volume", "Platelet_distribution_width", "High_light_scatter_reticulocyte_count",
               "UKB_460K.blood_LYMPHOCYTE_COUNT", "Monocyte_count", "UKB_460K.blood_MONOCYTE_COUNT", "UKB_460K.blood_RBC_DISTRIB_WIDTH")


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

  cis_trans_coloc <- fread(variables["lbf_based_cis_trans_coloc"]) %>%
    mutate(upstream_effect = T) %>%
    distinct(trans_gene, trans_locus, signal2, upstream_effect)

  gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
    filter(type == "gene") %>% select(gene_id, gene_name, seqid, start, end, strand) %>%
    mutate(TSS = if_else(strand == "+", start, end)) %>%
    filter(gene_id %in% available_genes)

  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type)) %>%
    select(variant, variant_index, phenotype, locus_str, SusieRss_CS, type, typeCisTrans, pip) %>%
    filter(!is.na(pip), phenotype %in% gene_reference$gene_id)

  # Process input
  variant_groups <- fread(variables["variant_coloc_classes"]) %>%
    separate_wider_delim("gene_id", delim="__", names=c("phenotype", "locus_str", "u_type"), too_many="drop", cols_remove = F) %>%
    mutate(SusieRss_CS = as.numeric(str_extract(cs_id, "\\d+")),
           locus_str = str_replace(locus_str, "_", ":"))

  independent_variant_set_annot_filtered <- left_join(independent_variant_set_annot_filtered, variant_groups)

  independent_variant_set_annot_filtered_with_cis_effects <-
    left_join(independent_variant_set_annot_filtered, variant_groups) %>%
    filter(typeCisTrans == "trans") %>%
    left_join(cis_trans_coloc, by = c("phenotype" = "trans_gene", "locus_str" = "trans_locus", "SusieRss_CS" = "signal2")) %>%
    mutate(upstream_effect = replace_na(upstream_effect, F))

  # Percentage of trans-eQTLs with upstream effects
  independent_variant_set_annot_filtered_with_cis_effects %>% count(upstream_effect) %>%
    mutate(p = n / sum(n) * 100)

  # Count number of eQTLs for which we identify GWAS colocalizations
  coloc_glob <- variables["coloc_cleaned"]

  colocalizations <- bind_rows(mapply(fread, Sys.glob(coloc_glob), SIMPLIFY = F)) %>%
    filter(eQTL_gene %in% gene_reference$gene_id) %>%
    group_by(eQTL_gene, signal1, eQTL_locus) %>%
    summarise(gwas_coloc = T,
           any_non_btc_coloc = any(!(GWAS %in% cc_traits)),
           all_non_btc_coloc = all(!(GWAS %in% cc_traits)),
           any_btc_coloc = any(GWAS %in% cc_traits))

  # Check enrichments of GWAS colocalization among trans-eQTLs that do have cis-effects vs those that do not.
  # ---
  # First perform enrichments on an eQTL level

  independent_variant_set_annot_filtered_with_cis_effects_flagged <- left_join(
    independent_variant_set_annot_filtered_with_cis_effects, colocalizations,
    by = c("phenotype" = "eQTL_gene", "SusieRss_CS" = "signal1", "locus_str" = "eQTL_locus")) %>%
    mutate(gwas_coloc = replace_na(gwas_coloc, F),
           any_btc_coloc = replace_na(any_btc_coloc, F),
           any_non_btc_coloc = replace_na(any_non_btc_coloc, F),
           all_non_btc_coloc_vs_no_coloc = case_when(
             all_non_btc_coloc ~ TRUE,
             is.na(all_non_btc_coloc) ~ FALSE,
             any_btc_coloc ~ NA
           ),
           all_non_btc_coloc = replace_na(all_non_btc_coloc, F))

  # Plot GWAS coloc enrichment
  ggplot(independent_variant_set_annot_filtered_with_cis_effects_flagged,
         aes(x=upstream_effect, fill=gwas_coloc)) +
    geom_bar(position = "fill")

  # Plot GWAS coloc enrichment excluding any blood cell type coloc
  ggplot(independent_variant_set_annot_filtered_with_cis_effects_flagged %>% filter(!is.na(all_non_btc_coloc_vs_no_coloc)),
         aes(x=upstream_effect, fill=all_non_btc_coloc_vs_no_coloc)) +
    geom_bar(position = "fill")

  # Test GWAS coloc enrichment
  fisher_out <- fisher.test(table(independent_variant_set_annot_filtered_with_cis_effects_flagged[,c("gwas_coloc", "upstream_effect")]))
  print(fisher_out)

  # Test GWAS coloc enrichment excluding any blood cell type coloc
  fisher_out <- fisher.test(table(independent_variant_set_annot_filtered_with_cis_effects_flagged[,c("all_non_btc_coloc_vs_no_coloc", "upstream_effect")]))
  print(fisher_out)

  # ---
  # Second, perform enrichments on an eQTL-variant level (this is presented in the manuscript)

  # Summarise eQTLs to independent loci
  independent_cs_set_annot_filtered_with_cis_effects_flagged <- independent_variant_set_annot_filtered_with_cis_effects_flagged %>%
    group_by(CS_group) %>%
    summarise(
      any_non_btc_coloc = any(any_non_btc_coloc),
      any_signal_btc_coloc = any(any_btc_coloc),
      all_non_btc_coloc = case_when(
        any(any_btc_coloc) ~ F,
        any(all_non_btc_coloc) ~ T,
        !any(gwas_coloc) ~ F
      ),
      all2_non_btc_coloc_vs_no_coloc = case_when(
        any_signal_btc_coloc ~ NA,
        any(all_non_btc_coloc) ~ T,
        !any(gwas_coloc) ~ F
      ),
      any_gwas_coloc = any(gwas_coloc),
      any_upstream_effect = any(upstream_effect)
    )

  # Calculate the number of trans-eQTL variants with cis-effects
  independent_cs_set_annot_filtered_with_cis_effects_flagged %>%
    count(any_upstream_effect) %>%
    ungroup() %>%
    mutate(p = n / sum(n))

  # Plot GWAS coloc enrichment
  ggplot(independent_cs_set_annot_filtered_with_cis_effects_flagged,
         aes(x=any_upstream_effect, fill=any_gwas_coloc)) +
    geom_bar(position = "fill")

  # Plot GWAS coloc enrichment excluding any blood cell type coloc
  ggplot(independent_cs_set_annot_filtered_with_cis_effects_flagged %>% filter(!is.na(all2_non_btc_coloc_vs_no_coloc)),
         aes(x=any_upstream_effect, fill=all2_non_btc_coloc_vs_no_coloc)) +
    geom_bar(position = "fill")

  # Test GWAS coloc enrichment
  fisher_out <- fisher.test(table(independent_cs_set_annot_filtered_with_cis_effects_flagged[,c("any_gwas_coloc", "any_upstream_effect")]))
  print(fisher_out)
  print(fisher_out$p.value)

  # Test GWAS coloc enrichment excluding any blood cell type coloc
  fisher_out <- fisher.test(table(independent_cs_set_annot_filtered_with_cis_effects_flagged[,c("all2_non_btc_coloc_vs_no_coloc", "any_upstream_effect")]))
  print(fisher_out)
  print(fisher_out$p.value)

  # ---
  # Now, plot the in-degree, out-degree distributions
  cis_trans_relationships <- fread(variables["cis_trans_relationships"])
  cis_trans_network <- graph_from_data_frame(cis_trans_relationships)

  fit_power_law(degree(cis_trans_network, mode="out"), p.value=T)
  fit_power_law(degree(cis_trans_network), p.value=T)

  deg_out <- degree(cis_trans_network, mode="out")

  # Compute P(k)
  deg_out_table <- table(deg_out)
  deg_out_dist <- data.frame(
    k = as.numeric(names(deg_out_table)),
    Pk = as.numeric(deg_out_table) / length(deg_out)
  )

  deg_in <- degree(cis_trans_network, mode="in")

  # Compute P(k)
  deg_in_table <- table(deg_in)
  deg_in_dist <- data.frame(
    k = as.numeric(names(deg_in_table)),
    Pk = as.numeric(deg_in_table) / length(deg_in)
  )

  deg_in_out <- bind_rows(list("in"=deg_in_dist, "out"=deg_out_dist), .id="mode")

  ggplot(deg_in_out, aes(x = k, y = Pk)) +
    geom_point(size = 1, alpha = 0.5, aes(shape=mode),  color="grey30") +
    scale_x_log10(limits=c(1, 10000)) +
    scale_y_log10(limits=c(1e-4, NA)) +
    scale_shape_manual(values=c(17, 16)) +
    labs(
      x = "Degree k (log scale)",
      y = "P(k) (log scale)"
    ) +
    theme(aspect.ratio = 1.0, panel.grid.major = element_line(color="grey90"))

  ggsave("degree_dist.pdf", width=2, height = 2)

  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}