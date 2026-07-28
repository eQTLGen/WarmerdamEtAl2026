#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)
library(igraph)
library(poweRlaw)
library(broom)
library(arrow)


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

  available_genes <- fread(variables["available_genes"], header=F)$V1
  variant_reference <- read_parquet(variables["variant_reference"])

  cis_trans_coloc <- fread(variables["lbf_based_cis_trans_coloc"]) %>%
    mutate(upstream_effect = T)

  gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
    filter(type == "gene") %>% select(gene_id, gene_name, seqid, start, end, strand) %>%
    mutate(TSS = if_else(strand == "+", start, end)) %>%
    filter(gene_id %in% available_genes)

  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type)) %>%
    select(variant, variant_index, phenotype, locus_str, SusieRss_CS, type, typeCisTrans, pip) %>%
    filter(!is.na(SusieRss_CS), phenotype %in% gene_reference$gene_id)

  cis_trans_coloc %>% inner_join(independent_variant_set_annot_filtered %>% select(phenotype, variant, type),  by = c("cis_gene" = "phenotype", "hit1" = "variant")) %>%
    filter(type == "cis")

  cis_trans_coloc %>% inner_join(independent_variant_set_annot_filtered %>% select(phenotype, variant, type),  by = c("cis_gene" = "phenotype", "hit1" = "variant")) %>%
    filter(type == "cis") %>% distinct(cis_gene, hit1)

  cis_trans_coloc_filtered <- cis_trans_coloc %>%
    inner_join(independent_variant_set_annot_filtered %>% select(phenotype, variant, type), by = c("cis_gene" = "phenotype", "hit1" = "variant")) %>%
    filter(type == "cis")

  cis_trans_coloc_filtered %>%
    distinct(trans_gene, signal2, trans_locus)

  # Process input
  variant_groups <- fread(variables["independent_credible_sets"]) %>%
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
    summarise(
      gwas_coloc = T,
      any_non_btc_coloc = any(!(GWAS %in% cc_traits)),
      all_non_btc_coloc = all(!(GWAS %in% cc_traits)),
      any_btc_coloc = any(GWAS %in% cc_traits))

  cis_trans_coloc_filtered %>% distinct(cis_locus, cis_gene, signal1) %>%
    left_join(colocalizations, by = c("cis_gene" = "eQTL_gene", "signal1" = "signal1", "cis_locus" = "eQTL_locus")) %>%
    summarise(n = n(),
              across(any_of(c("gwas_coloc", "any_non_btc_coloc", "all_non_btc_coloc", "any_btc_coloc")), function(x) sum(replace_na(x, FALSE))))

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

  ggplot(independent_variant_set_annot_filtered_with_cis_effects_flagged,
         aes(x=upstream_effect, fill=gwas_coloc)) +
    geom_bar(position = "fill")

  ggplot(independent_variant_set_annot_filtered_with_cis_effects_flagged,
         aes(x=upstream_effect, fill=any_non_btc_coloc)) +
    geom_bar(position = "fill")

  ggplot(independent_variant_set_annot_filtered_with_cis_effects_flagged,
         aes(x=upstream_effect, fill=all_non_btc_coloc)) +
    geom_bar(position = "fill")

  ggplot(independent_variant_set_annot_filtered_with_cis_effects_flagged %>% filter(!is.na(all_non_btc_coloc_vs_no_coloc)),
         aes(x=upstream_effect, fill=all_non_btc_coloc_vs_no_coloc)) +
    geom_bar(position = "fill")

  fisher_out <- fisher.test(table(independent_variant_set_annot_filtered_with_cis_effects_flagged[,c("gwas_coloc", "upstream_effect")]))
  print(fisher_out)

  fisher_out <- fisher.test(table(independent_variant_set_annot_filtered_with_cis_effects_flagged[,c("any_non_btc_coloc", "upstream_effect")]))
  print(fisher_out)

  fisher_out <- fisher.test(table(independent_variant_set_annot_filtered_with_cis_effects_flagged[,c("all_non_btc_coloc", "upstream_effect")]))
  print(fisher_out)

  fisher_out <- fisher.test(table(independent_variant_set_annot_filtered_with_cis_effects_flagged[,c("all_non_btc_coloc_vs_no_coloc", "upstream_effect")]))
  print(fisher_out)

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

  independent_cs_set_annot_filtered_with_cis_effects_flagged %>%
    count(any_upstream_effect) %>%
    ungroup() %>%
    mutate(p = n / sum(n))

  ggplot(independent_cs_set_annot_filtered_with_cis_effects_flagged,
         aes(x=any_upstream_effect, fill=any_gwas_coloc)) +
    geom_bar(position = "fill")

  ggplot(independent_cs_set_annot_filtered_with_cis_effects_flagged,
         aes(x=any_upstream_effect, fill=any_non_btc_coloc)) +
    geom_bar(position = "fill")

  ggplot(independent_cs_set_annot_filtered_with_cis_effects_flagged,
         aes(x=any_upstream_effect, fill=all_non_btc_coloc)) +
    geom_bar(position = "fill")

  ggplot(independent_cs_set_annot_filtered_with_cis_effects_flagged,
         aes(x=any_upstream_effect, fill=any_signal_btc_coloc)) +
    geom_bar(position = "fill")

  ggplot(independent_cs_set_annot_filtered_with_cis_effects_flagged %>% filter(!is.na(all2_non_btc_coloc_vs_no_coloc)),
         aes(x=any_upstream_effect, fill=all2_non_btc_coloc_vs_no_coloc)) +
    geom_bar(position = "fill")

  fisher_out <- fisher.test(table(independent_cs_set_annot_filtered_with_cis_effects_flagged[,c("any_gwas_coloc", "any_upstream_effect")]))
  print(fisher_out)
  print(fisher_out$p.value)

  fisher_out <- fisher.test(table(independent_cs_set_annot_filtered_with_cis_effects_flagged[,c("any_non_btc_coloc", "any_upstream_effect")]))
  print(fisher_out)
  print(fisher_out$p.value)

  fisher_out <- fisher.test(table(independent_cs_set_annot_filtered_with_cis_effects_flagged[,c("all_non_btc_coloc", "any_upstream_effect")]))
  print(fisher_out)

  fisher_out <- fisher.test(table(independent_cs_set_annot_filtered_with_cis_effects_flagged[,c("all2_non_btc_coloc_vs_no_coloc", "any_upstream_effect")]))
  print(fisher_out)
  print(fisher_out$p.value)

  cis_trans_relationships <- fread(variables["cis_trans_relationships"])
  mismatched_pairs <- fread(variables["cis_trans_pairs_mismatching"], fill=T)

  all_pairs_with_mismatches <- bind_rows("concordant" = cis_trans_relationships,"disconcordant" = mismatched_pairs, .id = "concordance") %>%
    filter(NrCredibleSetsAffectingBothCisAndTrans > 1) %>%
    mutate(p_null = 2*(0.5^NrCredibleSetsAffectingBothCisAndTrans))
  all_pairs_with_mismatches %>% group_by(NrCredibleSetsAffectingBothCisAndTrans, concordance) %>%
    tally() %>%
    group_by(NrCredibleSetsAffectingBothCisAndTrans) %>%
    summarise(p = n[concordance == "concordant"] / sum(n), sum(n))

  n <- all_pairs_with_mismatches$NrCredibleSetsAffectingBothCisAndTrans
  S_obs <- sum(all_pairs_with_mismatches$concordance == "concordant")

  f <- function(e) {
    sum((1 - e)^n + e^n) - S_obs
  }

  proportion_of_false_effect_directions <- uniroot(f, c(0, 0.5))$root
  # This is 3.6.

  # Observed count
  S_obs <- sum(all_pairs_with_mismatches$concordance == "concordant")

  # Expected count under H0
  E_S <- sum(all_pairs_with_mismatches$p_null)

  # Variance under H0
  Var_S <- sum(all_pairs_with_mismatches$p_null * (1 - all_pairs_with_mismatches$p_null))

  # Z statistic
  z <- (S_obs - E_S) / sqrt(Var_S)

  # One-sided p-value:
  # H1 = more same-color groups than expected
  p_value <- pnorm(z, lower.tail = F)

  cat("Observed same-color groups:", S_obs, "\n")
  cat("Expected under chance:", round(E_S, 3), "\n")
  cat("Z statistic:", round(z, 3), "\n")
  cat("One-sided p-value:", p_value, "\n")

  if (p_value < 0.05) {
    cat("\nEvidence for MORE same-color groups than expected by chance.\n")
  } else {
    cat("\nNo significant evidence for excess same-color groups.\n")
  }


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

  # Work
  deg <- degree(g, mode = "all")
  deg <- deg[deg > 0]

  deg_df <- data.frame(k = deg)
  ccdf <- function(x) { ecdf(x)(x) }

  ggplot(deg_df, aes(k)) +
    stat_ecdf(geom = "step") +
    scale_x_log10() +
    scale_y_log10(trans = "reverse") +  # makes it CCDF
    labs(x = "Degree k", y = "P(K ≥ k)", title = "Degree CCDF") +
    theme_minimal()


  m <- displ$new(deg)
  est <- estimate_xmin(m)
  m$setXmin(est)

  bs <- bootstrap_p(m, no_of_sims = 2000, threads=8)
  bs$p

  m_lognorm <- dislnorm$new(deg)
  est_lognorm <- estimate_xmin(m_lognorm)
  m_lognorm$setXmin(est_lognorm)

  comp <- compare_distributions(m, m_lognorm)
  comp

  m_exp <- disexp$new(deg)
  compare_distributions(m, m_exp)

  cis_trans_relationships <- fread(variables["cis_trans_relationships"])

  potential_cis_genes <- gene_reference %>%
    mutate(upstream_gene = gene_id %in% cis_trans_relationships$CisGene)

  gwas_types <- read_excel(variables["eqtlgen_gwass"]) %>%
    select(trait_identifier, trait_category) %>% deframe()

  coloc_cleaned <- fread(variables["coloc_cleaned"]) %>%
    filter(eQTL_type == "cis") %>%
    mutate(blood_trait_coloc = gwas_types[GWAS] == "Cell type composition") %>%
    group_by(eQTL_gene) %>%
    summarise(blood_trait_coloc = any(blood_trait_coloc))



  variant_classes <- fread(variables["variant_coloc_classes"]) %>% select(phenotype, variant_index, variant_type, type) %>%
    mutate(coloc = variant_type != "Not GWAS colocalising") %>%
    filter(type == "cis") %>%
    group_by(phenotype) %>%
    summarise(blood_trait_coloc = any(variant_type %in% c("Exclusively blood-cell", "Some blood-cell")))

  potential_cis_genes_with_coloc_annot <- potential_cis_genes %>%
    left_join(coloc_cleaned, by = c("gene_id" = "phenotype"))

  test_input <- potential_cis_genes_with_coloc_annot %>%
    mutate(blood_trait_coloc = replace_na(blood_trait_coloc, F)) %>%
    group_by(upstream_gene, blood_trait_coloc) %>% tally()

  test_input

  phyper(test_input$n[test_input$upstream_gene & test_input$blood_trait_coloc],
         sum(test_input$n[test_input$upstream_gene], na.rm=T),
         sum(test_input$n[!test_input$upstream_gene], na.rm=T),
         sum(test_input$n[test_input$blood_trait_coloc], na.rm=T),
         lower.tail = FALSE, log.p = FALSE)

  fisher.test(matrix(
    c(test_input$n[test_input$upstream_gene & test_input$blood_trait_coloc],
      test_input$n[test_input$upstream_gene & !test_input$blood_trait_coloc],
      test_input$n[!test_input$upstream_gene & test_input$blood_trait_coloc],
      test_input$n[!test_input$upstream_gene & !test_input$blood_trait_coloc]), nrow=2),
  )$p.value


  toppgene_mat <- fread(variables["toppgene_results"]) %>%
    as_tibble(.name_repair = c("unique")) %>%
    filter(Category=="GO: Molecular Function") %>%
    pivot_longer(cols = -c(Category, ID, Name, `Verbose ID`, pValue), names_to = "Gene", values_to = "Intersection")

  go_mat <- toppgene_mat %>%
    dplyr::select(ID, Name, pValue, Gene, Intersection) %>%
    pivot_wider(
      names_from = Gene,
      values_from = Intersection,
      values_fill = 0
    )

  # Keep metadata separately
  go_meta <- go_mat %>%
    dplyr::select(ID, Name, pValue)

  # Matrix for overlap calculations
  go_binary <- go_mat %>%
    dplyr::select(-ID, -Name, -pValue) %>%
    as.matrix()

  rownames(go_binary) <- go_mat$ID
  # Process output

  jaccard <- function(x, y) {
    sum(x & y) / sum(x | y)
  }

  n <- nrow(go_binary)

  jac <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(rownames(go_binary),
                    rownames(go_binary))
  )

  for(i in seq_len(n)){
    for(j in seq_len(n)){
      jac[i,j] <- jaccard(go_binary[i,], go_binary[j,])
    }
  }

  hc <- hclust(
    as.dist(1 - jac),
    method = "average"
  )

  mis_jac <- maximum_independent_set(jac > 0.1, rownames(jac))

  deg_out <- degree(cis_trans_network, mode = "out")

  cis_stats <- data.frame(
    gene = names(deg_out),
    out_degree = deg_out
  )

  cis_stats$pagerank <-
    page_rank(cis_trans_network, directed = TRUE)$vector[cis_stats$gene]

  cis_stats$coreness <-
    coreness(cis_trans_network,
             mode="all")[cis_stats$gene]

  cis_stats$harmonic <-
    harmonic_centrality(cis_trans_network,
                        mode="out")[cis_stats$gene]

  cis_stats$betweenness <-
    betweenness(cis_trans_network, normalized=T, directed = TRUE)[cis_stats$gene]

  cis_stats$closeness <-
    closeness(cis_trans_network, mode = "out", normalized = T)[cis_stats$gene]

  hgnc_symbols <- fread(variables["hgnc_symbols"]) %>% dplyr::select(ensembl_gene_id, symbol) %>% deframe()

  # Compare classes:
  toppfun <- fread("/Users/c.a.warmerdam/Downloads/2a2f6ffb-cad1-4836-9558-ac3bd0928a31.txt") %>%
    separate_longer_delim(`Hit in Query List`, delim=",")

  cis_stats_per_enrichment_gene_set <- cis_stats %>%
    mutate(hgnc_symbol = hgnc_symbols[gene]) %>%
    left_join(toppfun, by=c("hgnc_symbol" = "Hit in Query List"))

  ggplot(cis_stats_per_enrichment_gene_set %>% filter(ID %in% mis_jac$independent_vertices), aes(y=Name, x=out_degree)) +
    geom_violin() + geom_boxplot(width=0.1) + scale_x_log10()

  ggplot(cis_stats_per_enrichment_gene_set %>% filter(ID %in% mis_jac$independent_vertices), aes(y=Name, x=betweenness)) +
    geom_violin() + geom_boxplot(width=0.1) + scale_x_log10()

  ggplot(cis_stats_per_enrichment_gene_set %>% filter(ID %in% mis_jac$independent_vertices), aes(y=Name, x=closeness)) +
    geom_violin() + geom_boxplot(width=0.1)

  dat <- cis_stats_per_enrichment_gene_set %>%
    filter(ID %in% mis_jac$independent_vertices) %>%
    filter(!is.na(Name))

  library(rstatix)

  kruskal_test(pagerank ~ Name, data = dat)

  kruskal_test(harmonic ~ Name, data = dat)

  kruskal_test(coreness ~ Name, data = dat)

  kruskal_test(out_degree ~ Name, data = dat)

  kruskal_test(betweenness ~ Name, data = dat)

  kruskal_test(closeness ~ Name, data = dat)

  #This gives a single p-value:
  # "Network centrality differs across functional classes."

  # 2. Pairwise comparisons
  # If the global test is significant:
  test_class <- function(class_name, metric, dat) {

    dat2 <- dat %>%
      mutate(
        group = ifelse(Name == class_name,
                       class_name,
                       "Other")
      )

    wt <- wilcox.test(
      reformulate("group", metric),
      data = dat2
    )

    class_values <- dat2 %>%
      filter(group == class_name) %>%
      pull(all_of(metric))

    other_values <- dat2 %>%
      filter(group == "Other") %>%
      pull(all_of(metric))

    tibble(
      class = class_name,
      metric = metric,
      n_class = length(class_values),
      n_other = length(other_values),
      median_class = median(class_values, na.rm=TRUE),
      median_other = median(other_values, na.rm=TRUE),
      fold_change = median_class / median_other,
      p = wt$p.value
    )
  }

  classes <- unique(dat$Name)
  metrics <- c(
    "out_degree",
    "harmonic",
    "coreness"
  )

  test_results <- expand_grid(
    class = classes,
    metric = metrics
  ) %>%
    mutate(
      result = pmap(
        list(class, metric),
        ~test_class(..1, ..2, dat)
      )
    ) %>%
    select(result) %>%
    unnest(result) %>%
    group_by(metric) %>%
    mutate(
      p_adj = p.adjust(p, method="BH")
    ) %>%
    ungroup()

  library(ggbeeswarm)

  plot_metric <- function(metric, dat, test_results){

    # global test
    kw <- kruskal_test(
      reformulate("Name", metric),
      data = dat
    )

    # order classes by median
    dat_plot <- dat %>%
      mutate(
        Name = fct_reorder(
          Name,
          .data[[metric]],
          median,
          na.rm=TRUE
        )
      )


    # labels with n
    labels <- dat_plot %>%
      count(Name) %>%
      mutate(
        label = paste0(
          Name,
          "\n(n=",
          n,
          ")"
        )
      )


    # significant class tests
    annotations <- test_results %>%
      filter(
        metric == !!metric,
        p < 0.05
      ) %>%
      mutate(
        Name = factor(
          class,
          levels = levels(dat_plot$Name)
        ),
        label = paste0(
          "P=",
          format.pval(
            p,
            digits=2,
            eps=0.001
          )
        ),
        x = max(
          dat_plot[[!!metric]],
          na.rm=TRUE
        ) * 1
      )


    p <- ggplot(
      dat_plot,
      aes(
        x=.data[[metric]],
        y=Name
      )
    ) +
      geom_violin(
        fill="grey90",
        scale="width"
      ) +
      geom_quasirandom(
        width=0.15,
        alpha=0.5,
        size=1
      ) +
      geom_boxplot(
        width=0.12,
        outlier.shape=NA
      ) +
      scale_y_discrete(
        labels=setNames(
          labels$label,
          labels$Name
        )
      ) +
      labs(
        x=metric,
        y=NULL,
        title=paste0(
          metric,
          " differs between functional classes"
        ),
        subtitle=paste0(
          "Kruskal-Wallis p=",
          signif(kw$p,3)
        )
      ) +
      theme_classic()


    # log scale for degree-like metrics
    if(metric %in% c("out_degree")){
      p <- p +
        scale_x_log10()
    }


    # Add significant labels
    if(nrow(annotations)>0){

      p <- p +
        geom_text(
          data=annotations,
          aes(
            x=x,
            y=Name,
            label=label
          ),
          inherit.aes=FALSE,
          hjust=1,
          size=3
        ) +
        coord_cartesian(
          clip="off"
        )
    }

    p
  }

  plot_metric(
    "out_degree",
    dat,
    test_results
  )

  plot_metric(
    "harmonic",
    dat,
    test_results
  )

  plot_metric(
    "coreness",
    dat,
    test_results
  )
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}