#!/usr/bin/env Rscript
#` Code for generating forest plots in supplements

# Load libraries
library(tidyverse)
library(arrow)
library(data.table)
library(igraph)

# Declare constants
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


# Gene combinations
gene_gene_combinations <- c("EEPD1" = "SREBF1",
                            "EEPD1" = "ABCA1",
                            "EEPD1" = "SQLE",
                            "SLC39A10" = "SLC30A1",
                            "SLC39A10" = "MT1X",
                            "SLC39A10" = "MT2A",
                            "SLC39A10" = "MT1F",
                            "STAT2" = "OAS2",
                            "STAT2" = "OAS3",
                            "STAT2" = "IFI27",
                            "STAT2" = "IFI44",
                            "STAT2" = "SIGLEC1",
                            "STAT2" = "USP18"
)


compute_paths_to_focal_gene <- function(g, path_table, focal_gwas_id, focal_gene_id) {
  path_table_2 <- path_table %>% filter(gene == focal_gene_id, gwas_id == focal_gwas_id)
  association_names <- unique(path_table_2 %>% pull(association))

  # Edge direction vector (assumed to be numeric: 1 or -1)
  edge_directions <- E(g)$Direction

  results <- list()

  for (s in association_names) {
    # Get all simple paths up to 2 steps
    paths <- all_simple_paths(g, from = s, to = focal_gene_id, cutoff = 2)

    for (p in paths) {
      node_ids <- as_ids(p)
      print(node_ids)

      # Get first edge
      eid1 <- get_edge_ids(g, node_ids[1:2])
      path_sign <- edge_directions[eid1]

      if (length(node_ids) > 2) {
        eid2 <- get_edge_ids(g, node_ids[2:3])
        path_sign <- path_sign * edge_directions[eid2]
      }

      edge_attrs <- E(g)[eid1]

      tibbleResult <- tibble(
        source = node_ids[1],
        target = node_ids[length(node_ids)],
        path_length = length(node_ids) - 1,
        path_sign = path_sign,
        eqtl_variant = edge_attrs$hit1,
        first_edge_type = edge_attrs$edgeType,
        first_edge_cistrans = edge_attrs$typeCisTrans,
        first_edge_id = eid1,
        intermediate_node = ifelse(length(node_ids) == 3, node_ids[2], NA_character_)
      )
      results[[length(results) + 1]] <- tibbleResult
      print(tibbleResult)
    }
  }

  print(results)

  bound_paths <- bind_rows(results) %>%
    filter(!(path_length == 2 & first_edge_cistrans == "trans"))
  return(bound_paths)
}

forest_plot <- function(per_cohort_data, meta_analysed_data, hgnc_map, variant_reference, master_table, focal_variant, genes, plot_legend = F) {
  per_cohort_results <- per_cohort_data %>% filter(variant_index == !!focal_variant) %>% collect() %>%
    filter(phenotype %in% genes)
  meta_results <- meta_analysed_data %>% filter(variant_index == !!focal_variant) %>% collect() %>%
    filter(phenotype %in% genes) %>% mutate(cohort = "meta")

  print("Extracting results done!")
  print("cisGene")
  print(hgnc_map[genes])
  cisGene <- hgnc_map[genes[1]]
  print(cisGene)

  meta_and_per_cohort_bound <- bind_rows(per_cohort_results, meta_results) %>%
    inner_join(variant_reference) %>%
    mutate(symbol = hgnc_map[phenotype],
           symbol = factor(symbol, levels=unique(c(cisGene, symbol)), ordered=T)) %>%
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
      legend.position = if_else(plot_legend, "right", "none")
    )

  return(p)
}

generate_forest_plots <- function(
  per_cohort_data, meta_analysed_data, hgnc_map, variant_reference, master_table, example_variant_gene_combinations) {
  return(lapply(example_variant_gene_combinations %>% group_by(variant_index) %>% group_split(), function(group) {
    print(group)
    focal_variant <- group$variant_index[1]
    forest_plot(per_cohort_data, meta_analysed_data, hgnc_map,
                variant_reference = variant_reference, master_table = master_table,
                focal_variant = focal_variant, genes = unlist(group$genes[1]), plot_legend = focal_variant == example_variant_gene_combinations$variant_index[1])
  }))
}


# Declare function definitions
extract_gene_eqtls <- function(per_cohort_data, meta_data, gene, variant_reference) {
  cohort_data <- per_cohort_data %>% filter(phenotype == !!gene) %>% collect()
  meta_data <- meta_analysed %>% filter(phenotype == !!gene) %>% select(variant_index, beta, standard_error, phenotype) %>% mutate(cohort="meta")
  meta_and_per_cohort_bound <- bind_rows(cohort_data, meta_data) %>% inner_join(variant_reference)
}

get_abca1_network <- function() {

  focal_gene <- "ABCA1"
  focal_gwas_id <- "Low_Density_Lipoprotein_Cholesterol_levels"
  focal_gene_id <- hgnc_genes[symbol == focal_gene, ]$ensembl_gene_id

  prefix <- "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/network_analyses/revision/coloc_run_20260602_AllEqtls_UnLd05_AllGwas"
  g <- readRDS(paste0(prefix, "_graph.rds"))

  association_nodes_df <- tibble(igraph::as_data_frame(g, what = "vertices")[V(g)$nodeType == "association", ]) %>%
    select(name, gwas_id, index_variant) %>%
    mutate(variant = case_when(
      !is.na(index_variant) ~ index_variant,
      TRUE ~ str_replace(name, "[^:]+\\:", "")
    ))

  path_table_filtered <- readRDS(paste0(prefix, "_path_table_processed.rds")) %>%
    inner_join(association_nodes_df, by = c("association" = "name"))

  paths_to_focal_gene <- compute_paths_to_focal_gene(g, path_table_filtered, focal_gwas_id, focal_gene_id)
  abca1_network <- paths_to_focal_gene %>% distinct(source, target, intermediate_node, eqtl_variant) %>%
    inner_join(association_nodes_df, by = c("source" = "name")) %>%
    distinct(variant, eqtl_variant, intermediate_node, target)
}

# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  resources <- "resources/roberts_mbp.tsv"
  resources <- "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/users/umcg-rwarmerdam/resources.txt"
  variables <- deframe(read_delim(resources, delim="\t"))
  location <- str_split(variables["available_genes"], "freeze3")[[1]][1]

  hgnc_genes <- fread("/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/network_analyses/hgnc_complete_set.txt")
  hgnc_map <- deframe(hgnc_genes %>% select(ensembl_gene_id, symbol))

  cohort_info <- fread("/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/cohort_data/cohort_info.txt")
  master_table <- fread(variables["master_table"]) %>% select(cohort, cohort_new_name, dataset_n, dataset_platform, dataset_ancestry) %>%
    inner_join(cohort_info, by = c("cohort_new_name" = "cohort"))

  variant_reference <- read_parquet(variables["variant_reference"])

  example_gene_gene_pairs <- enframe(gene_gene_combinations)

  causal_gene_gene_relationships <- fread(variables["cis_trans_relationships"])
  example_gene_gene_relationships <- causal_gene_gene_relationships %>%
    mutate(CisGeneName = hgnc_map[CisGene],
           TransGeneName = hgnc_map[TransGene]) %>%
    semi_join(example_gene_gene_pairs, by = c("CisGeneName" = "name", "TransGeneName" = "value"))

  example_genes <- unique(c(example_gene_gene_relationships$CisGene, example_gene_gene_relationships$TransGene))

  example_variant_gene_combinations <- example_gene_gene_relationships %>%
    group_by(CisSNPs) %>%
    summarise(genes = list(unique(c(CisGene, TransGene)))) %>%
    inner_join(variant_reference, by = c("CisSNPs" = "variant"))

  abca1_network <- get_abca1_network()

  abca1_genes <- na.omit(unique(abca1_network$intermediate_node, abca1_network$target))

  abca1_variant_gene_combinations <- abca1_network %>%
    group_by(eqtl_variant) %>%
    summarise(genes = list(unique(c(intermediate_node, target)))) %>%
    inner_join(variant_reference, by = c("eqtl_variant" = "variant"))

  all_genes_to_sync <- unique(c(example_genes, abca1_genes))
  fwrite(
    tibble(folders_to_sync = sprintf("phenotype=%s", all_genes_to_sync)),
    "folders_to_sync.txt",
    col.names=F, row.names=F
  )

  # Process input
  variant_reference <- read_parquet(variables["variant_reference"])

  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"])

  per_cohort_data <- open_dataset("/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/eqtl_mapping/output/all_empirical_per_cohort_4GenPC20ExpPC_2025-06-12/eqtls/cohort")
  meta_analysed_data <- open_dataset("/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/eqtl_mapping/output/all_empirical_per_cohort_4GenPC20ExpPC_2025-06-12/eqtls/meta")

  example_gene_forest_plots <- generate_forest_plots(
    per_cohort_data, meta_analysed_data, hgnc_map,
    variant_reference, master_table,
    example_variant_gene_combinations)

  pdf("example_gene_forest_plots_20260722.pdf", onefile = TRUE, width=180/10/2.54, height = 180/10/2.54)
  for (i in seq(length(example_gene_forest_plots))) {
    print(example_gene_forest_plots[[i]])
  }
  dev.off()

  abca1_network_forest_plots <- generate_forest_plots(
    per_cohort_data, meta_analysed_data, hgnc_map,
    variant_reference, master_table,
    abca1_variant_gene_combinations)

  pdf("abca1_network_forest_plots_20260722.pdf", onefile = TRUE, width=180/10/2.54, height = 180/10/2.54)
  for (i in seq(length(abca1_network_forest_plots))) {
    print(abca1_network_forest_plots[[i]])
  }
  dev.off()

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}