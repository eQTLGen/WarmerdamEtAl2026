#!/usr/bin/env Rscript
#` This code summarises output of the gwas gene prioritization

# Load libraries
library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(emmeans)
library(broom)

# Declare constants
old <- theme_set(theme_classic(base_size = 8))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  panel.grid.major.y = element_line(color="grey70"),
  panel.grid.minor = element_blank(),
  text = element_text(family="Helvetica"),
  title = element_text(colour = "#595A5C", face = "plain")
)


# Declare function definitions

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

  with_edges <- colSums(adjacency_matrix[unconnected_vertices, active_nodes, drop=F] == T) > 0
  active_nodes[with_edges] <- F

  clusters <- list()

  # Loop while the matrix is not empty
  while (any(active_nodes)) {
    # Find the vertex with the least amount of edges
    least_connected <- names(base::which.min(colSums(adjacency_matrix[active_nodes, active_nodes, drop=F])))

    # Add the least connected vertex to the independent set
    independent_vertices <- c(independent_vertices, least_connected)

    # Remove all vertices that have an edge with the least connected vertex
    with_edges <- adjacency_matrix[least_connected,][active_nodes, drop=F] == T

    active_nodes[names(with_edges)[with_edges]] <- F

    clusters[[least_connected]] <- names(active_nodes[names(with_edges)[with_edges]])
  }

  return(list(independent_vertices=independent_vertices, clusters=clusters))
}


prune_susie_output <- function(chromosome_eqtl_table, ld_func, r2_threshold) {

  if (nrow(chromosome_eqtl_table) == 1) {
    return(chromosome_eqtl_table)
  }

  # Get the LD matrix
  ld_matrix <- ld_func(chromosome_eqtl_table$chromosome[1], chromosome_eqtl_table$variant_index)
  ld_matrix <- ld_matrix[
    as.character(chromosome_eqtl_table$variant_index),
    as.character(chromosome_eqtl_table$variant_index)]

  if (!all(rownames(ld_matrix) == as.character(chromosome_eqtl_table$variant_index))) {
    stop("LD matrix row names must match variant IDs in variants table")
  }

  MIS_out <- maximum_independent_set((ld_matrix^2) > r2_threshold, rownames(ld_matrix))

  clusters <- enframe(MIS_out$clusters) %>%
    mutate(variant_index = as.numeric(name)) %>%
    select(variant_index, value)

  chromosome_eqtl_table_pruned <- chromosome_eqtl_table %>% filter(as.character(variant_index) %in% MIS_out$independent_vertices) %>%
    left_join(clusters, by = c("variant_index" = "variant_index")) %>%


    return(chromosome_eqtl_table_pruned)
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

  gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
    filter(type == "gene") %>% select(gene_id, gene_name, seqid, start, end, strand)

  eqtls <- fread(variables["independent_variant_set"]) %>%
    filter(!is.na(pip)) %>%
    select(variant_index, phenotype, chromosome, type)
  cis_trans_pairs <- fread(variables["cis_trans_relationships"]) %>%
    select(CisGene, TransGene)

  all_possible_paths <-
    eqtls %>% filter(type == "cis") %>% inner_join(cis_trans_pairs, by = c("phenotype" = "CisGene")) %>%
      select(peripheral = phenotype, phenotype = TransGene, variant_index, chromosome, type) %>%
      bind_rows(eqtls) %>%
      group_by(phenotype) %>%
      distinct(variant_index, .keep_all = T)

  all_possible_paths %>% filter(phenotype == "ENSG00000204388")

  eqtl_matrix <- as.matrix(fread(variables["z_mat"]), rownames = "variant")
  sign_threshold_nominal <- qnorm(0.05*2, lower.tail=F)

  idx <- cbind(as.character(all_possible_paths$variant_index), all_possible_paths$phenotype)
  all_possible_paths$eQTL_z <- eqtl_matrix[idx]
  all_possible_paths$nominal <- abs(all_possible_paths$eQTL_z) > sign_threshold_nominal

  precalculated_ld_path <- variables["finemapped_variant_ld"]
  if (file.exists(precalculated_ld_path)) {
    ld_matrices <- readRDS(precalculated_ld_path)
  }

  ld_func_mat <- function(chromosome, variant_indices) {
    chromosome_ld_matrix <- ld_matrices[[as.character(chromosome)]]
    return(chromosome_ld_matrix[as.character(variant_indices), as.character(variant_indices)])
  }

  output_prune_function <- bind_rows(mapply(function(.x) {
    prune_susie_output(.x, ld_func_mat, 0.5)
  }, all_possible_paths %>% filter(nominal) %>%
                                              group_by(chromosome, phenotype) %>% group_split(), SIMPLIFY=F))

  saveRDS(output_prune_function, "theoretical_paths.rds")

  gwas_selection <- readxl::read_xlsx(variables["eqtlgen_gwass"]) %>%
    select(trait_identifier, trait_name, trait_category, type) %>%
    mutate(
      trait_category_fine = case_when(
        trait_category %in% c("Cell type composition", "Immune mediated diseases") ~ trait_category,
        trait_category == "Metabolism & cardiometabolic disease" & type == "quant" ~ "Metabolism & cardiometabolic traits",
        trait_category == "Metabolism & cardiometabolic disease" & type == "cc" ~ "Metabolism & cardiometabolic diseases",
        trait_category == "Neuropsychiatric & cognitive traits" & type == "quant" ~ "Neuropsychiatric & cognitive traits",
        trait_category == "Neuropsychiatric & cognitive traits" & type == "cc" ~ "Neuropsychiatric & cognitive diseases",
        trait_category == "Other traits" & type == "cc" ~ "Other diseases",
        trait_category == "Other traits" & type == "quant" ~ "Other traits",
      ),
      trait_category_alt = case_when(
        trait_category == "Cell type composition" ~ trait_category,
        type == "cc" ~ "Diseases",
        type == "quant" ~ "Quantitative traits"
      ),
      trait_category_alt2 = case_when(
        trait_category == "Cell type composition" ~ trait_category,
        trait_category %in% c("Immune mediated diseases", "Metabolism & cardiometabolic disease") & type == "cc" ~ "Blood-related diseases",
        trait_category %in% c("Immune mediated diseases", "Metabolism & cardiometabolic disease") & type == "quant" ~ "Blood-related quantitative traits")
    )

  # Fill in!
  output_prefix <- "/Users/c.a.warmerdam/Documents/projects/2025/eQTLGen/freeze3/Interpretation/network_analyses/revision/coloc_run_20260602_AllEqtls_UnLd05_AllGwas"

  # Process input
  gene_prioritizations_all <- fread(sprintf("%s_prioritization_perm_results.rds", output_prefix)) %>%
    mutate(gene_score_p_adj_b = p.adjust(gene_score_p, method = "bonferroni"), gene_score_p_adj_bh = p.adjust(gene_score_p, method = "BH"),
           fdr_significant = gene_score_p_adj_bh < 0.05, Bonferroni_significant = gene_score_p_adj_b < 0.05) %>%
    left_join(gwas_selection, by=c("gwas_id" = "trait_identifier")) %>%
    select(-trait_category_alt2, -trait_category_alt, -trait_category_fine)

  fwrite(gene_prioritizations_all, paste0(output_prefix, "gene_prioritization_results_all.txt"), sep="\t", row.names=F, col.names=T)

  gene_prioritizations_all %>%
    filter(gene_score_p_adj_bh < 0.05, trait_category != "Cell type composition") %>%
    count(abs(gene_score) > 1) %>%
    mutate(n / sum(n))

  gene_prioritizations_all %>%
    filter(gene_score_p_adj_bh < 0.05, trait_category != "Cell type composition", abs(gene_score) <= 1)

  significant_results <- gene_prioritizations_all %>%
    filter(gene_score_p_adj_bh < 0.05, abs(gene_score) > 1, trait_category != "Cell type composition")

  significant_results %>%
    group_by(sign(gene_score)) %>%
    summarise(n = n()) %>%
    mutate(p = n / sum(n) * 100)

  significant_results %>%
    group_by(gwas_id) %>%
    summarise(n=n()) %>%
    summarise(n_studies = n(), sum = sum(n), mean = sum / n_studies, mean_all = sum / 70)
  # Perform method
  # Process output

  theoretical_paths <- readRDS(variables["theoretical_paths_network_analyses"])
  max_finemapped_eqtls_per_gene_ld_05 <- theoretical_paths %>%
    group_by(phenotype) %>%
    summarise(n_eqtls = n())

  gene_prioritizations_all_extended <- fread(sprintf("%s_prioritization_perm_results.rds", output_prefix)) %>%
    mutate(gene_score_p_adj_b = p.adjust(gene_score_p, method = "bonferroni"), gene_score_p_adj_bh = p.adjust(gene_score_p, method = "BH"),
           fdr_significant = gene_score_p_adj_bh < 0.05, Bonferroni_significant = gene_score_p_adj_b < 0.05)

  gene_prioritizations_all_extended_b <- expand.grid(gene = gene_reference$gene_id, gwas_id = unique(gene_prioritizations_all_extended$gwas_id)) %>%
    anti_join(gene_prioritizations_all_extended, by = c("gene", "gwas_id")) %>%
    bind_rows(gene_prioritizations_all_extended) %>%
    left_join(gwas_selection, by=c("gwas_id" = "trait_identifier")) %>%
    left_join(max_finemapped_eqtls_per_gene_ld_05, by = c("gene" = "phenotype"))

  repressing_cols <- grep("^repressing_", names(gene_prioritizations_all_extended_b))
  activating_cols <- grep("^activating_", names(gene_prioritizations_all_extended_b))

  gene_prioritizations_all_extended_2 <- gene_prioritizations_all_extended_b %>%
    mutate(
      gene_score = if_else(is.na(gene_score), 0, gene_score),
      sum_repressing = rowSums(across(all_of(repressing_cols)), na.rm = TRUE),
      sum_activating = rowSums(across(all_of(activating_cols)), na.rm = TRUE),
      sum_activating_or_repressing = pmax(sum_repressing, sum_activating),
      network_size = sum_repressing + sum_activating,
      gene_score_vs_network_size = abs(gene_score) / network_size,
      large_network = network_size > 1,
      opposing_direction = abs(network_size - abs(gene_score)) > 0,
      network_size_vs_eqtls = network_size / n_eqtls,
      opposing_directions = (network_size - abs(gene_score)) / 2
    )

  gene_prioritizations_all_extended_2 %>%
    filter(fdr_significant, trait_category != "Cell type composition") %>%
    ungroup() %>%
    filter(abs(gene_score) <= 1) %>%
    group_by(network_size) %>%
    summarise(n = n(),
              across(starts_with("sign_sum"), function(x) sum(x != 0, na.rm=T)),
    ) %>%
    mutate(p = n / sum(n))

  gene_prioritizations_all_extended_2 %>%
    filter(fdr_significant, trait_category != "Cell type composition") %>%
    group_by(trait_category) %>%
    summarise(mean(network_size), median(network_size), mean(network_size_vs_eqtls), mean(n_eqtls), n(), sum(network_size > 1), sum(n_eqtls > 1))

  gene_prioritizations_all_extended_2 %>%
    filter(fdr_significant, trait_category != "Cell type composition") %>%
    group_by(gene_score_abs = abs(gene_score), network_size) %>%
    summarise(mean(network_size), median(network_size), mean(network_size_vs_eqtls), mean(n_eqtls), median(n_eqtls), n(), sum(network_size > 1), sum(n_eqtls > 1)) %>%
    mutate(p = n / sum(n))


  trait_summary <- gene_prioritizations_all_extended_2 %>%
    group_by(trait_name, trait_category, trait_category_fine, trait_category_alt, trait_category_alt2, type) %>%
    summarise(
      n_opposing = sum(opposing_direction, na.rm=T),
      n_opposing_sig = sum(opposing_direction[fdr_significant], na.rm=T),
      n_large = sum(network_size > 1, na.rm=T),
      n_large_sig = sum(network_size[fdr_significant] > 1, na.rm=T),
      prop_large_sig = mean(network_size[fdr_significant] > 1, na.rm=T),
      prop_large = mean(network_size > 1),
    )


  comparisons <- tribble(
    ~group_var,            ~group1,                    ~group2,
    "trait_category_alt",  "Cell type composition",    "Diseases",
    "trait_category_alt",  "Cell type composition",    "Quantitative traits",
    "trait_category_alt2", "Cell type composition",    "Blood-related diseases",
    "trait_category_alt2", "Cell type composition",    "Blood-related quantitative traits"
  )

  outcomes <- c("prop_large", "prop_large_sig")

  results <- list()

  for (i in seq_len(nrow(comparisons))) {

    comp <- comparisons[i, ]

    for (outcome in outcomes) {

      x <- trait_summary[[outcome]][
        trait_summary[[comp$group_var]] == comp$group1
      ]

      y <- trait_summary[[outcome]][
        trait_summary[[comp$group_var]] == comp$group2
      ]

      test_res <- broom::tidy(
        wilcox.test(na.omit(x), na.omit(y))
      )

      results[[length(results) + 1]] <- bind_cols(
        tibble(outcome = outcome,
               grouping_variable = comp$group_var,
               group_1 = comp$group1,
               group_2 = comp$group2,
               median_group_1 = median(x, na.rm = T),
               median_group_2 = median(y, na.rm = T),
               mean_group_1 = mean(x, na.rm = T),
               mean_group_2 = mean(y, na.rm = T)),
        test_res
      )
    }
  }

  results <- bind_rows(results)

  ggplot(trait_summary, aes(trait_category_alt, prop_large)) + geom_violin() + geom_boxplot(width=0.1) + scale_y_log10() + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))
  ggplot(trait_summary, aes(trait_category_alt2, prop_large)) + geom_violin() + geom_boxplot(width=0.1) + scale_y_log10() + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))

  my_breaks <- c(1, 10, 100, 1000, 10000)

  p1 <- gene_prioritizations_all_extended_2 %>%
    filter(fdr_significant, trait_category != "Cell type composition") %>%
    mutate(gene_score_abs = factor(abs(gene_score), ordered=T)) %>%
    group_by(trait_category, network_size, gene_score_abs, opposing_directions) %>%
    tally() %>% filter(network_size > 0) %>%
    ggplot(aes(x=network_size, y=opposing_directions, fill=gene_score_abs)) +
    geom_tile() + geom_text(aes(label=n)) + facet_wrap(vars(trait_category), scales="free") +
    scale_fill_viridis_d() +
    scale_x_continuous(breaks=c(0,1,2,3,4,5,6)) +
    scale_y_continuous(breaks=c(0,1,2)) +
    ylab("N variants with an opposing direction") +
    xlab("Network size in number of linked variants")

  p2 <- gene_prioritizations_all_extended_2 %>%
    filter(fdr_significant, trait_category == "Cell type composition") %>%
    mutate(gene_score_abs = abs(gene_score)) %>%
    group_by(trait_category, network_size, gene_score_abs, opposing_directions) %>%
    tally() %>% filter(network_size > 0) %>%
    ggplot(aes(x=network_size, y=opposing_directions, fill=gene_score_abs)) +
    geom_tile() + geom_text(aes(label=n)) + facet_wrap(vars(trait_category), scales="free") +
    scale_fill_viridis_c() +
    ylab("N variants with an opposing direction") +
    xlab("Network size in number of linked variants")

  p_comp <- grid.arrange(p2, p1, ncol=1)

  ggsave("compiled.png", p_comp, dpi=300, height=180, width=180, units="mm")


  trait_summary2 <- gene_prioritizations_all_extended_2 %>%
    group_by(trait_name, trait_category, trait_category_fine, trait_category_alt, trait_category_alt2, type) %>%
    filter(abs(network_size) == 3) %>%
    summarise(
      average_ratio = mean(((network_size - abs(gene_score)) / 2) / network_size),
      n = n()
    )

  ggplot(trait_summary2, aes(x=trait_category_alt2, y=average_ratio)) + geom_violin() + geom_boxplot(width=0.1)

  kruskal.test(average_ratio ~ trait_category_alt, data = trait_summary2)

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}