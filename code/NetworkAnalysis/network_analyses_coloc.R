#!/usr/bin/env Rscript


# Load libraries
library(argparse)
library(data.table)
library(tidyverse)
library(dtplyr)
library(igraph)
library(extrafont)
library(arrow)
library(ggraph)
library(tidygraph)
library(GenomicRanges)
library(parallel)
library(rtracklayer)

# Declare constants
old <- theme_set(theme_classic(base_size = 8))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  panel.grid.minor = element_blank(),
  text = element_text(family="Helvetica"),
  title = element_text(colour = "#595A5C", face = "plain")
)

# Declare constants
gwas_ids_of_interest <- c(
  "PASS_Celiac", "PASS_IBD_deLange2017", "PASS_Rheumatoid_Arthritis", "PASS_Lupus_Wang2021",
  "PASS_Coronary_Artery_Disease", "PASS_CD_deLange2017", "PASS_Multiple_sclerosis2019",
  "PASS_UC_deLange2017", "PASS_AtrialFibrillation_Nielsen2018", "UKB_460K.disease_CARDIOVASCULAR", "UKB_460K.disease_ASTHMA_DIAGNOSED", "PASS_Type_1_Diabetes")

gwas_ids_of_interest <- c(
  "PASS_IBD_deLange2017", "UKB_460K.disease_CARDIOVASCULAR", "UKB_460K.disease_ASTHMA_DIAGNOSED", "PASS_Lupus_Wang2021")

gwas_labels <- c(
  PASS_Celiac = "Celiac Disease",
  PASS_IBD_deLange2017 = "Inflammatory Bowel Disease",
  PASS_Lupus = "Lupus",
  PASS_Lupus_Wang2021 = "SLE (Wang et al. 2021)",
  PASS_Rheumatoid_Arthritis = "Rheumatoid Arthritis",
  PASS_Coronary_Artery_Disease = "Coronary Artery Disease",
  PASS_CD_deLange2017 = "Crohn's Disease",
  PASS_Multiple_sclerosis2019 = "Multiple Sclerosis",
  PASS_UC_deLange2017 = "Ulcerative Colitis",
  PASS_AtrialFibrillation_Nielsen2018 = "Atrial Fibrillation",
  UKB_460K.disease_CARDIOVASCULAR = "Cardiovascular Disease",
  UKB_460K.disease_ASTHMA_DIAGNOSED = "Asthma",
  PASS_Type_1_Diabetes = "Type 1 Diabetes"
)

non_blood_cell_type_traits <- c("PASS_CD_deLange2017", "PASS_IBD_deLange2017",
                                "PASS_Lupus_Wang2021", "PASS_Multiple_sclerosis2019", "PASS_Primary_biliary_cirrhosis",
                                "PASS_Rheumatoid_Arthritis", "PASS_UC_deLange2017", "UKB_460K.disease_AID_ALL",
                                "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED",
                                "UKB_460K.disease_ASTHMA_DIAGNOSED", "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP",
                                "UKB_460K.disease_RESPIRATORY_ENT", "PASS_ADHD_Demontis2023",
                                "PASS_Alzheimers_Bellenguez2022", "PASS_BIP_Mullins2021", "PASS_DrinksPerWeek_Liu2019",
                                "PASS_GeneralRiskTolerance_KarlssonLinner2019", "PASS_Insomnia_Watanabe2022",
                                "PASS_Intelligence_SavageJansen2018", "PASS_MDD_Als2023", "PASS_ReactionTime_Davies2018",
                                "PASS_SWB", "PASS_Schizophrenia_Trubetskoy2022", "PASS_SleepDuration_Dashti2019",
                                "PASS_VerbalNumericReasoning_Davies2018", "PASS_Worry_Nagel2018",
                                "UKB_460K.body_BMIz", "UKB_460K.cov_EDU_COLLEGE", "UKB_460K.cov_EDU_YEARS",
                                "UKB_460K.cov_SMOKING_STATUS", "UKB_460K.mental_NEUROTICISM",
                                "UKB_460K.other_MORNINGPERSON", "UKB_460K.repro_NumberChildrenEverBorn_Pooled",
                                "PASS_AtrialFibrillation_Nielsen2018", "PASS_Coronary_Artery_Disease",
                                "UKB_460K.bp_DIASTOLICadjMEDz", "UKB_460K.bp_SYSTOLICadjMEDz",
                                "UKB_460K.disease_CARDIOVASCULAR", "UKB_460K.disease_HYPERTENSION_DIAGNOSED",
                                "PASS_FastingGlucose_Manning", "UKB_460K.biochemistry_AlanineAminotransferase",
                                "UKB_460K.biochemistry_AlkalinePhosphatase", "UKB_460K.biochemistry_Cholesterol",
                                "UKB_460K.biochemistry_Glucose", "UKB_460K.biochemistry_HDLcholesterol",
                                "UKB_460K.biochemistry_HbA1c", "UKB_460K.biochemistry_LDLdirect",
                                "UKB_460K.biochemistry_SHBG", "UKB_460K.biochemistry_Testosterone_Male",
                                "UKB_460K.biochemistry_TotalBilirubin", "UKB_460K.biochemistry_TotalProtein",
                                "UKB_460K.biochemistry_Triglycerides", "PASS_Type_1_Diabetes",
                                "PASS_Type_2_Diabetes", "UKB_460K.bmd_HEEL_TSCOREz", "UKB_460K.body_BALDING1",
                                "UKB_460K.body_HEIGHTz", "UKB_460K.body_WHRadjBMIz", "UKB_460K.cancer_BREAST",
                                "UKB_460K.impedance_BASAL_METABOLIC_RATEz", "UKB_460K.lung_FEV1FVCzSMOKE",
                                "UKB_460K.lung_FVCzSMOKE", "UKB_460K.pigment_HAIR", "UKB_460K.repro_MENARCHE_AGE",
                                "UKB_460K.repro_MENOPAUSE_AGE", "High_Density_Lipoprotein_Cholesterol_levels",
                                "Low_Density_Lipoprotein_Cholesterol_levels", "Non-HDL_Cholesterol_levels",
                                "Total_Cholesterol_Measurement", "Triglyceride_Levels")

# Declare function definitions
save_pdf <- function(p, gene, gwas) {
  file.pathResult <- file.path(".", sprintf("gene_gene_subgraph.%s_%s.%s.pdf", gwas, gene, format(Sys.Date(), "%Y%m%d")))
  pdf(
    file.pathResult, width = 12, height=12)

  print(p)

  dev.off()
}


calculate_gene_scores <- function(disease_gene_paths) {

  # Summarise sign_sum
  gene_scores <- disease_gene_paths[
    , .(sign_sum = sum(sign),
        activating = sum(sign>0),
        repressing = sum(sign<0),
        associations = paste(eqtl_variant, sign, collapse=",", sep=":")),
      by = .(gene, dist, first_edge_cis_trans)
  ]

  # Pivot to wide format
  gene_scores_wide <- dcast(
    gene_scores,
    gene ~ first_edge_cis_trans + dist,
    value.var = c("sign_sum", "activating", "repressing", "associations"),
    fill = 0
  )

  # Compute gene_score as rowSums of all numeric columns except "gene"
  gene_scores_wide[
    , gene_score := rowSums(.SD),
      .SDcols = names(gene_scores_wide)[str_detect(names(gene_scores_wide), "^sign_sum_")]
  ]

  # Return result
  return(gene_scores_wide)
}

random_variant_permutation_test <- function(g, focal_gwas_id, path_table_filtered, gwas_dist_mat) {
  association_nodes <- which(V(g)$nodeType == "association" & V(g)$gwas_id == focal_gwas_id)
  focal_associations <- V(g)$name[association_nodes]
  cat(sprintf("Found %s focal associations", length(focal_associations)))

  association_nodes_df <- tibble(as_data_frame(g, what = "vertices")[V(g)$nodeType == "association", ]) %>%
    select(name, gwas_id)

  path_table_filtered_with_gwas_id <- inner_join(path_table_filtered, association_nodes_df, by = c("association" = "name"))

  if (length(focal_associations) < 1) {
    return(NULL)
  }

  uncorrelated_to_focal_gwas <- names(which(!gwas_dist_mat[focal_gwas_id,]))

  uncorrelated_gwas_set_including_focal_gwas <- maximum_independent_set(
    gwas_dist_mat[uncorrelated_to_focal_gwas, uncorrelated_to_focal_gwas], uncorrelated_to_focal_gwas)

  gwas_subset <- c(uncorrelated_gwas_set_including_focal_gwas$independent_vertices, focal_gwas_id)
  cat(sprintf("Found %s independent GWASs including focal one", length(gwas_subset)))

  # Filter valid paths
  disease_gene_paths <- path_table_filtered_with_gwas_id[
    association %in% focal_associations
  ]

  actual_gene_scores <- calculate_gene_scores(disease_gene_paths)

  # Step 1: Get n_distinct_associations per group
  n_associations <- disease_gene_paths[
    , .(n_distinct_associations = uniqueN(association)),
      by = .(first_edge_cis_trans, dist)
  ]

  # Step 2: Join with valid_paths and filter by gwas_subset
  permutation_sampling_data <- merge(
    n_associations,
    path_table_filtered_with_gwas_id,
    by = c("first_edge_cis_trans", "dist")
  )[gwas_id %in% gwas_subset]

  calculate_permuted_scores <- function(iteration) {
    cat(sprintf("\rProgress: %d", iteration))

    # Step 1: Sample associations WITH replacement per group
    sampled_associations <- permutation_sampling_data[
      , {
          # Unique associations in the group
          unique_assocs <- unique(association)

          # Sample associations with replacement
          sampled <- sample(
            unique_assocs,
            size = n_distinct_associations[1],
            replace = TRUE
          )

          .(association = sampled)
        },
        by = .(first_edge_cis_trans, dist)
    ]

    # Step 2: Join back to get all corresponding rows (including duplicates)
    sampled_paths <- sampled_associations[
      permutation_sampling_data,
      on = .(first_edge_cis_trans, dist, association),
      nomatch = 0
    ]

    return(list(
      gene_scores=calculate_gene_scores(sampled_paths)))
  }

  message("")

  permutations <- lapply(1:1000, calculate_permuted_scores)

  message("Done!")

  permuted_gene_scores <- bind_rows(lapply(permutations, function(x) x$gene_scores), .id = "iteration")

  # Define the score columns of interest
  score_columns <- c("gene_score")

  # Compute mean and sd for each gene and each score column in the background
  background_stats <- permuted_gene_scores %>%
    group_by(gene) %>%
    summarise(
      zeros = 1000-n(),
      across(all_of(score_columns),
             list(mean = ~mean(c(.x, rep(0, zeros)), na.rm = TRUE),
                  sd = ~sd(c(.x, rep(0, zeros)), na.rm = TRUE)),
             .names = "{.col}_{.fn}"))

  # Join the background stats with the actual gene scores
  gene_scores_with_stats <- actual_gene_scores %>%
    dplyr::select("gene", any_of(c(score_columns)), ends_with("cis_1"), ends_with("cis_2"), ends_with("trans_1")) %>%
    left_join(background_stats, by = "gene") %>%
    mutate(across(all_of(score_columns),
                  ~ if_else(is.na(get(paste0(cur_column(), "_mean"))), Inf,
                            (. - get(paste0(cur_column(), "_mean"))) / get(paste0(cur_column(), "_sd"))),
                  .names = "{.col}_z")) %>%
    mutate(across(all_of(score_columns),
                  ~ (2*pnorm(abs(get(paste0(cur_column(), "_z"))), lower.tail = F)),
                  .names = "{.col}_p"))

  return(gene_scores_with_stats)
}


tToCor <- function(t, df){
  if(length(t) != length(df)){
    return(0L)  #Error condition
  }
  r <- t/sqrt(df+t^2)
  return(r)
}


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

plot_focused_network2 <- function(g, focal_gwas_id, focal_gene_id, paths_to_focal_gene, hgnc_gene_names, bands_map) {

  # Get all nodes from paths
  relevant_nodes <- na.omit(c(paths_to_focal_gene$source, paths_to_focal_gene$intermediate_node, focal_gene_id))

  # Define subgraph with all relevant nodes and edges
  g_sub <- induced_subgraph(g, vids = relevant_nodes)

  gene_score <- paths_to_focal_gene %>%
    group_by(source) %>%
    slice_min(path_length, with_ties=F) %>%
    ungroup() %>%
    summarise(sum_of_path_signs = sum(path_sign))

  # Node attributes: color and size by node type
  node_fill_map <- c(
    "focal" = "white",
    "gene" = "white",
    "association" = "#B3B3B3"
  )

  node_color_map <- c(
    "focal" = "black",
    "gene" = "black",
    "association" = "#666666"
  )

  node_size_map <- c(
    "focal" = 4,
    "gene" = 2,
    "association" = 1
  )

  # Edge line types and colors
  edge_linetype_map <- c(
    "genepairs" = "solid",
    "cis" = "solid",
    "trans" = "longdash"
  )

  edge_color_map <- c(
    "1" = "#5AB847",    # Activating
    "-1" = "#E41E25"    # Repressing
  )

  # Node type levels (for ordering)
  node_type_levels <- c("focal", "gene", "association")

  # Node label mappings
  label_map <- name_map[focal_gene_id]
  title_text <- sprintf("%s - %s", gwas_labels[focal_gwas_id], label_map)

  # ------------------------------------------------------------------------------
  # Construct tidygraph object with node and edge attributes
  # ------------------------------------------------------------------------------

  g_tbl <- as_tbl_graph(g_sub) %>%
    mutate(
      # Assign node types
      node_type = case_when(
        name == focal_gene_id ~ "focal",
        nodeType == "gene" ~ "gene",
        nodeType == "association" ~ "association"
      ),

      # Assign bands
      bands = case_when(
        nodeType == "gene" ~ bands_map[name],
        nodeType == "association" ~ bands_map[variant]
      ),

      # Assign node labels
      label = case_when(
        name %in% names(hgnc_gene_names) ~ sprintf("%s", hgnc_gene_names[name]),
        node_type %in% c("association") ~ bands_map[variant],
        node_type == "proxy" ~ variant
      ),

      # Assign colors and font styles
      color = case_when(
        label == label_map ~ node_color_map["focal"],
        node_type == "gene" ~ node_color_map["gene"],
        node_type == "association" ~ node_color_map["association"]
      ),
      fontface = case_when(
        node_type %in% c("focal", "gene") ~ "italic",
        TRUE ~ "plain"
      ),

      # Factorize node type for layout ordering
      node_type_fctr = ordered(node_type, levels = node_type_levels)
    ) %>%
    arrange(desc(node_type_fctr)) %>%
    activate(edges) %>%
    mutate(
      # Assign edge types and properties
      edge_type = case_when(
        edgeType == "genelink" ~ "genepairs",
        typeCisTrans == "cis" ~ "cis",
        typeCisTrans == "trans" ~ "trans"
      ),
      edge_sign = Direction,
      edge_weight = if_else(edge_type == "ld", 1, 10),
      arrow_type = "closed"
    )

  # ------------------------------------------------------------------------------
  # Generate layout with edge weights
  # ------------------------------------------------------------------------------

  layout <- create_layout(g_tbl, layout = "kk", weights = edge_weight)

  # ------------------------------------------------------------------------------
  # Plot graph using ggraph
  # ------------------------------------------------------------------------------

  p <- ggraph(layout) +

    # Positive regulatory edges
    geom_edge_link(
      data = get_edges("short")(layout) %>% filter(edge_sign == 1, edge_type != "ld"),
      aes(color = as.character(Direction), linetype = edge_type),
      linewidth = 0.32,
      show.legend = TRUE,
      arrow = arrow(length = unit(1.2, 'mm')),
      start_cap = circle(2, 'mm'),
      end_cap = circle(2, 'mm')
    ) +

    # Negative regulatory edges
    geom_edge_link(
      data = get_edges("short")(layout) %>% filter(edge_sign == -1, edge_type != "ld"),
      aes(color = as.character(Direction), linetype = edge_type),
      linewidth = 0.32,
      show.legend = TRUE,
      arrow = arrow(length = unit(0.72, 'mm'), angle = 90),
      start_cap = circle(2, 'mm'),
      end_cap = circle(2, 'mm')
    ) +

    # Nodes: points and text
    geom_node_point(aes(color = node_type, size = node_type, fill=node_type), show.legend = TRUE, shape=21) +

    # Chromosome band labels for variants and ld_clusters
    geom_node_text(
      data = layout %>% filter(node_type %in% c("gene", "focal")),
      aes(label=label, fontface = fontface),
      repel = TRUE,
      size = 2.2,
      force = 20,
      vjust = -1.2,
      color = "black",
      show.legend = FALSE
    ) +

    # Scales
    scale_edge_linetype_manual(values = edge_linetype_map, na.value = "solid") +
    scale_edge_color_manual(values = edge_color_map, labels = c("1" = "Activating", "-1" = "Repressing")) +
    scale_color_manual(values = node_color_map, labels = c(
      "focal" = "Prioritized gene",
      "gene" = "Gene",
      "association" = "GWAS Association"
    )) +
    scale_fill_manual(values = node_fill_map, labels = c(
      "focal" = "Prioritized gene",
      "gene" = "Gene",
      "association" = "GWAS Association"
    )) +
    scale_size_manual(values = node_size_map) +

    # Theme and title
    theme_void() +
    theme(legend.position = "none") +
    ggtitle(title_text)

  # Chromosome band labels for variants and ld_clusters
  p <- p +
    geom_node_text(
      data = layout %>% filter(node_type %in% c("association")),
      aes(label = bands),
      repel = TRUE,
      size = 2.2,
      force = 20,
      vjust = -1.5,
      color = "black",
      show.legend = FALSE
    )

  # Save the plot
  ggsave(paste0(name_map[focal_gene_id], "_", focal_gwas_id, "_2025-07-09.pdf"), p, width = 6/2.54, height = 5/2.54)
}


get_distance_matrix <- function(g, cached_path) {
  if (file.exists(cached_path)) {
    return(readRDS(cached_path))
  }

  association_nodes <- which(V(g)$nodeType == "association")
  gene_nodes <- which(V(g)$nodeType == "gene")
  gene_names <- V(g)$name[gene_nodes]
  association_names <- V(g)$name[association_nodes]

  # Distance matrix: rows = SNPs, cols = genes
  dist_matrix <- distances(g, v = association_names, to = gene_names, mode = "out")

  saveRDS(dist_matrix, cached_path)

  return(dist_matrix)
}


calculate_sign_matrix <- function(g, dist_matrix, gene_nodes, distance_threshold) {

  # 1. Identify SNP and gene nodes
  association_nodes <- which(V(g)$nodeType == "association")
  if (is.null(gene_nodes)) {
    gene_nodes <- which(V(g)$nodeType == "gene")
  }

  gene_names <- V(g)$name[gene_nodes]
  association_names <- V(g)$name[association_nodes]

  # 2. Edge direction vector (assumed to be numeric: 1 or -1)
  edge_directions <- E(g)$Direction
  edge_type_cis_trans <- E(g)$typeCisTrans

  combs_that_pass_dist_threshold <- (dist_matrix <= distance_threshold)

  # 6. Prepare for parallel gene-wise computation
  compute_gene_column <- function(i) {
    gene_id <- gene_nodes[i]
    snp_indices <- which(combs_that_pass_dist_threshold[,i])

    if (length(snp_indices) == 0) return(NULL)

    snp_names_subset <- association_names[snp_indices]

    paths <- shortest_paths(g, from = gene_id, to = snp_names_subset, mode = "in", output = "epath")

    cat(sprintf("\rProgress: %d", i))

    out <- vector("list", length(snp_indices))

    for (jf in seq_along(snp_indices)) {
      path <- paths$epath[[jf]]
      if (length(path) == 0) next

      j <- unname(snp_indices[jf])

      sign_val <- prod(edge_directions[path])
      first_edge_cis_trans <- edge_type_cis_trans[path[length(path)]]

      first_edge_gene <- ends(g, rev(path)[1])[2] # From the node ids of the first edge, get the second node
      # This needs to be reversed because we use mode 'in' in the shortest paths thing

      out[[jf]] <- c(j = j, i = i, sign = sign_val, first_edge_cis_trans = first_edge_cis_trans, first_edge_gene = first_edge_gene)
    }

    return(out)
  }

  gene_results <- mclapply(seq_along(gene_names), compute_gene_column, mc.cores=16)
  saveRDS(gene_results, "path_calculations_max_2_edges.rds")

  gene_results_as_dt <- lapply(list_flatten(gene_results), function(mat) {
    if (!is.null(mat)) {
      dt <- as.data.table(t(mat))
    } else {
      dt <- data.table()
    }
    return(dt)
  })

  path_data_bound <- rbindlist(gene_results_as_dt)

  path_data_bound$i <- as.numeric(path_data_bound$i)
  path_data_bound$j <- as.numeric(path_data_bound$j)

  path_data_bound$gene <- gene_names[path_data_bound$i]
  path_data_bound$association <- association_names[path_data_bound$j]

  path_data_bound <- path_data_bound %>%
    rowwise() %>%
    mutate(dist = dist_matrix[j, i],
           sign = as.integer(sign))

  return(path_data_bound)
}

# Get sign matrix
get_path_table <- function(g, dist_matrix, gene_nodes, cached_path) {
  if (!is.na(cached_path) && cached_path != F && file.exists(cached_path)) {
    return(readRDS(cached_path))
  }
  path_table <- calculate_sign_matrix(g, dist_matrix, gene_nodes, 2)

  saveRDS(path_table, cached_path)
  return(path_table)
}


cs_clusters_per_gwas <- function(ld_matrices, gwas_eqtl_coloc_filtered_annot, gwas_ids = NULL, r2_threshold = 0.1) {
  if (is.null(gwas_ids)) {
    gwas_ids <- unique(gwas_eqtl_coloc_filtered_annot$GWAS)
  }

  variant_clusters_as_tib_all <- bind_rows(mapply(function(focal_gwas_id) {
    gwas_variant_indices <- gwas_eqtl_coloc_filtered_annot %>% filter(GWAS == focal_gwas_id) %>% pull(eqtl_variant_index)

    variant_clusters <- list_flatten(lapply(ld_matrices, function(ld_mat) {
      chr_gwas_variant_indices <- as.character(intersect(as.integer(rownames(ld_mat)), gwas_variant_indices))
      if (length(chr_gwas_variant_indices) <= 1) {
        return(list())
      }
      clusters <- maximum_independent_set(
        (ld_mat[chr_gwas_variant_indices, chr_gwas_variant_indices]^2) > r2_threshold,
        nodes=chr_gwas_variant_indices)$clusters

      return(Filter(function(x) length(x) > 1, clusters))
    }))

    variant_clusters_as_tib <- tibble(
      name = names(variant_clusters),
      values = lapply(variant_clusters, as.integer)
    )

    return(variant_clusters_as_tib)
  }, gwas_ids, SIMPLIFY=F, USE.NAMES = T), .id = "gwas_id")

  variant_clusters_as_tib_all <- variant_clusters_as_tib_all %>%
    separate_wider_delim(name, delim="_", names=c("chromosome", "index_variant_index")) %>%
    mutate(index_variant_index = as.integer(index_variant_index))

  return(variant_clusters_as_tib_all)
}


# Build network
build_network <- function(gwas_eqtl_coloc_filtered_annot, cis_trans, variant_annotations, ld_matrices, r2_threshold = 0.5) {
  variant_map <- deframe(variant_annotations %>% select(variant_index, variant))

  variant_clusters_as_tib_all <- cs_clusters_per_gwas(ld_matrices, gwas_eqtl_coloc_filtered_annot, r2_threshold = r2_threshold)
  # Create table of GWAS association (basis of nodes and edges from gwas associations)

  ld_clusters <- variant_clusters_as_tib_all %>%
    select(-chromosome) %>%
    inner_join(variant_annotations, by = c("index_variant_index" = "variant_index")) %>%
    group_by(gwas_id, chr_bands) %>%
    mutate(ld_cluster = paste(chr_bands, row_number(), sep="_ld-cluster_")) %>%
    ungroup() %>%
    select(gwas_id, ld_cluster, index_variant_index, values)

  # Create table of GWAS association (basis of nodes and edges from gwas associations)
  # Per GWAS, per ld thing, rename the hit variants
  gwas_associations <- ld_clusters %>%
    unnest_longer(col=values, values_to="eqtl_variant_index") %>%
    right_join(gwas_eqtl_coloc_filtered_annot, by = c("gwas_id"="GWAS", "eqtl_variant_index"="eqtl_variant_index")) %>%
    select(-cs95_PIP, -cs95_variants) %>%
    mutate(
      eqtl_variant = variant_map[as.character(eqtl_variant_index)],
      index_variant = variant_map[as.character(index_variant_index)],
      association_type = if_else(is.na(ld_cluster), "variant", "ld_cluster"),
      name = if_else(is.na(ld_cluster), paste0(gwas_id, ":", eqtl_variant), paste0(gwas_id, "_", ld_cluster))) %>%
    group_by(gwas_id, eQTL_gene, name) %>%
    slice_max(abs(eQTL_lead_eQTL_Z), with_ties = F, n=1)

  # Create GWAS associations-to-gene eQTL nodes
  # Here we filter to only those association-gene eqtls for which the GWAS normalized effect sizes have the same direction of effect
  gwas_gene_edges <- gwas_associations %>%
    mutate(
      from = name,
      to = eQTL_gene,
      eqtl_variant = eqtl_variant,
      eqtl_variant_index = eqtl_variant_index,
      edgeType = "eqtl",
      Direction_gwas = sign(eQTL_lead_GWAS_Z),
      sign_eqtl = sign(eQTL_lead_eQTL_Z),
      Direction = Direction_gwas * sign_eqtl,
      typeCisTrans = eQTL_type,
      edgeType = "eqtl",
      risk_allele_bulk_z = Direction_gwas * eQTL_lead_eQTL_Z,
      risk_allele_sc_z = Direction_gwas * sc_z,
      across(starts_with("sc_z_score_"), function(z) z * Direction_gwas)
    ) %>%
    group_by(from, to) %>%
    filter(n_distinct(Direction) == 1) %>%
    slice_head()

  # Direction: is the risk allele increasing expression of the eGene (+1), or decreasing the expression of the eGene (-1)

  # create GWAS association nodes
  gwas_association_nodes <- gwas_associations %>%
    semi_join(gwas_gene_edges, by = c("name"="from")) %>%
    group_by(name) %>%
    slice_head(n=1) %>%
    mutate(nodeType = "association") %>%
    select(name, nodeType, index_variant_index, index_variant, gwas_id)

  # Create gene-gene edges
  gene_gene_edges <- cis_trans %>% select(from = CisGene, to = TransGene, Direction) %>%
    mutate(edgeType = "genelink")

  gene_nodes <- bind_rows(
    tibble(name = c(gene_gene_edges$from, gene_gene_edges$to), nodeType = "gene"),
    tibble(name = gwas_gene_edges$to, nodeType = "gene")) %>%
    distinct()

  all_edges <- bind_rows(
    gene_gene_edges,
    gwas_gene_edges,
  )

  all_nodes <- bind_rows(
    gene_nodes,
    gwas_association_nodes,
  )

  fwrite(all_edges,
         sprintf("cis-trans_network_all_edges_coloc.txt.gz"),
         sep="\t", row.names=F, col.names=T, quote=F, na="NA")

  g <- graph_from_data_frame(all_edges, vertices = all_nodes, directed = TRUE)
}


filter_path_table <- function(path_table, g, eqtl_matrix) {

  # Per path, identify what the variant is that links the association to the first_edge_gene

  eqtl_edges <- as_data_frame(g, what = "edges")[E(g)[edgeType == "eqtl"], ] %>%
    select(association = from, first_edge_gene = to, eqtl_variant_index, eqtl_variant, Direction_gwas)

  # Path table is prefiltered to a max distance of 2.

  # Filter paths on those that are direct,
  # or those for which the first edge is a cis-eQTL and that is no longer than 2
  paths_category_filtered_with_gwas_beta <- path_table %>%
    as.data.table() %>%
    mutate(path_id = row_number(),
           sign = as.numeric(sign)) %>%
    filter(dist == 1 | first_edge_cis_trans == "cis") %>%
    select(path_id, sign, first_edge_gene, first_edge_cis_trans, gene, association, dist) %>%
    inner_join(eqtl_edges, by=c("association", "first_edge_gene"))

  # For these paths, require that the source SNP is associated to the focal gene with at least nominal significance
  sign_threshold_nominal <- qnorm(0.05*2, lower.tail=F)

  # Per path, get the eQTL gene
  # Per path, get the SNP from the LD cluster that is linked to the intermediary gene

  # Ensure variant_index is character (matches rownames of eqtl_matrix)
  variant_index <- as.character(paths_category_filtered_with_gwas_beta$eqtl_variant_index)
  gene <- paths_category_filtered_with_gwas_beta$gene

  # Get row and column indices
  row_idx <- match(variant_index, rownames(eqtl_matrix))
  col_idx <- match(gene, colnames(eqtl_matrix))

  # Use matrix subsetting to extract values
  direct_z_scores <- eqtl_matrix[cbind(row_idx, col_idx)]

  # Add results back into data frame
  filtered_paths_with_eqtl_z_score <- paths_category_filtered_with_gwas_beta %>%
    mutate(
      direct_z_score = direct_z_scores,
      nominal = abs(direct_z_score) > sign_threshold_nominal,
      sign_concordant = Direction_gwas * sign(direct_z_score) == sign)

  filtering_summary <- filtered_paths_with_eqtl_z_score %>%
    count(nominal, sign_concordant, dist, first_edge_cis_trans)
  print(tibble(filtering_summary), n=Inf)

  paths_category_filtered_and_nominal_significance_filtered <- filtered_paths_with_eqtl_z_score %>%
    filter(sign_concordant, nominal)

  return(paths_category_filtered_and_nominal_significance_filtered)
}


merge_overlapping <- function(lst) {
  n <- length(lst)
  adj <- matrix(FALSE, n, n)

  # build adjacency matrix: edges if overlap
  for (i in 1:n) {
    for (j in (i+1):n) {
      if (overlaps(lst[[i]], lst[[j]])) {
        adj[i, j] <- adj[j, i] <- TRUE
      }
    }
  }

  # find connected components
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  comps <- igraph::components(g)$membership

  # merge within components
  res <- lapply(split(seq_along(lst), comps), function(idx) {
    Reduce(union, lst[idx])
  })

  res
}


filter_eqtl_coloc_table <- function(gwas_eqtl_coloc, sc_eqtl_files, force_single_cell_replication) {

  # First, find eQTLs credible sets that colocalize to multiple loci of the same GWAS
  # gwas_eqtl_coloc_dedup <- deduplicate_gwas_loci(gwas_eqtl_coloc)

  names(sc_eqtl_files) <- paste0("sc_z_score_", str_extract(sc_eqtl_files, "\\.(\\w+)\\.incPerm\\.tsv$", group =1))
  sc_replication <- bind_rows(mapply(fread, sc_eqtl_files, SIMPLIFY=F, USE.NAMES=T), .id="celltype") %>%
    mutate(z_score = -1 * z_score) %>%
    pivot_wider(id_cols = c(alt_snp_id, alt_feature_id),values_from = "z_score", names_from = "celltype")

  # For these paths, require that the source SNP is associated to the focal gene with at least nominal significance
  sign_threshold_nominal <- qnorm(0.05*2, lower.tail=F)

  gwas_eqtl_coloc_with_sc <- gwas_eqtl_coloc %>%
    left_join(sc_replication, by = c("hit1" = "alt_snp_id", "eQTL_gene" = "alt_feature_id")) %>%
    mutate(
      any_ct_tested = if_any(starts_with("sc_z_score_"), ~ !is.na(.x)),
      # turn the 6 cols into a matrix
      any_ct_nominal_concordant = if_any(
        starts_with("sc_z_score_"),
        ~ replace_na((abs(.x) > sign_threshold_nominal) & (sign(.x) == sign(eQTL_lead_eQTL_Z)), F))
    )

  sc_summarised <- gwas_eqtl_coloc_with_sc %>%
    filter(any_ct_nominal_concordant) %>%
    # go long on the sc z-scores
    pivot_longer(
      cols = starts_with("sc_z_score_"),
      names_pattern = "(sc_z_score)_(.+)",
      names_to = c(".value", "cell_type")
    ) %>%
    # mark concordance
    mutate(
      sc_z = as.numeric(sc_z_score),
      is_concordant = !is.na(sc_z) &
        (abs(sc_z) > sign_threshold_nominal) &
        (sign(sc_z) == sign(eQTL_lead_eQTL_Z))
    ) %>%
    filter(is_concordant) %>%
    group_by(hit1, eQTL_gene) %>%
    # pick the most significant concordant cell type
    slice_max(order_by = abs(sc_z), n = 1, with_ties = FALSE) %>%
    select(
      sc_z, cell_type, hit1, eQTL_gene) %>%
    ungroup()

  gwas_eqtl_coloc_with_sc_annot <- gwas_eqtl_coloc_with_sc %>%
    left_join(sc_summarised, by = c("hit1", "eQTL_gene"))

  gwas_eqtl_coloc_with_sc_annot %>%
    count(any_ct_nominal_concordant)

  if (force_single_cell_replication) {

    return(gwas_eqtl_coloc_with_sc_annot %>% filter(any_ct_nominal_concordant))
  }

  return(gwas_eqtl_coloc_with_sc_annot)
}


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

get_gwas_distance_matrix <- function(gwas_dataset_path, uncorrelated_variant_path, cache_path = "pairwise_gwas_correlations.rds") {

  if (!file.exists(cache_path)) {

    if (is.null(gwas_ds)) {
      gwas_ds <- open_dataset(gwas_dataset_path)
    }

    uncorrelated_variants <- fread(uncorrelated_variant_path, header=F, col.names=c("variant")) %>%
      inner_join(variant_annotations)

    uncorrelated_gwas_variants <- gwas_ds %>%
      filter(variant_index %in% uncorrelated_variants$variant_index) %>%
      select(gwas_id, beta, se, variant_index, N) %>%
      collect()

    pearson_r_matrix <- uncorrelated_gwas_variants %>%
      mutate(pearson_rho = tToCor(beta/se, N - 1)) %>%
      group_by(gwas_id, variant_index) %>%
      slice_min(beta/se) %>%
      pivot_wider(id_cols="variant_index", names_from="gwas_id", values_from = "pearson_rho") %>%
      as.data.table() %>% as.matrix(rownames="variant_index")

    pairwise_gwas_correlations <- cor(pearson_r_matrix, use="pairwise.complete.obs")
    saveRDS(pairwise_gwas_correlations, cache_path)

  } else {
    pairwise_gwas_correlations <- readRDS(cache_path)
    pairwise_gwas_correlations <- pairwise_gwas_correlations[rownames(pairwise_gwas_correlations) != "PASS_Lupus",
                                                             colnames(pairwise_gwas_correlations) != "PASS_Lupus"]

  }
  gwas_dist_mat <- pairwise_gwas_correlations > 0.2
  return(gwas_dist_mat)
}


# Helper to load config TSV into a named list
load_config <- function(path) {
  deframe(read_delim(path, delim = "\t", col_names = c("key", "value")))
}


# Utility function: safely get config value, or NULL if missing
get_or_null <- function(config, key) {
  if (!is.null(config) && key %in% names(config)) {
    return(config[[key]])
  } else {
    return(NULL)
  }
}



# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  print(argv)

  # Phase 1: Get config path only
  parser0 <- ArgumentParser(add_help = FALSE)
  parser0$add_argument("--config", help = "Path to TSV config file", required = TRUE)
  args0 <- parser0$parse_known_args(argv)
  config_path <- args0[[1]]$config
  if (!file.exists(config_path)) stop("Config file not found: ", config_path)
  config <- load_config(config_path)

  print(config)

  # Phase 2: Build full parser with config values as defaults
  parser <- ArgumentParser(description = "Run gene prioritization pipeline")
  parser$add_argument("--config", help = "Path to config file")  # keep for completeness

  # Core inputs
  parser$add_argument("--hgnc_file", default = get_or_null(config, "hgnc_file"))
  parser$add_argument("--cytoband_file", default = get_or_null(config, "cytoband_file"))
  parser$add_argument("--variant_annotations", default = get_or_null(config, "variant_annotations"))
  parser$add_argument("--gene_reference", default = get_or_null(config, "gene_reference"))
  parser$add_argument("--cis_trans_file", default = get_or_null(config, "cis_trans_file"))
  parser$add_argument("--independent_variant_set", default = get_or_null(config, "independent_variant_set"))
  parser$add_argument("--ld_matrices", default = get_or_null(config, "ld_matrices"))
  parser$add_argument("--gwas_dataset", default = get_or_null(config, "gwas_dataset"))
  parser$add_argument("--uncorrelated_variant_path", default = get_or_null(config, "uncorrelated_variant_path"))
  parser$add_argument("--gwas_eqtl_coloc", default = get_or_null(config, "gwas_eqtl_coloc"))
  parser$add_argument("--eqtl_matrix", default = get_or_null(config, "eqtl_matrix"))
  parser$add_argument("--sc_replication", default = get_or_null(config, "sc_replication_path"))
  parser$add_argument("--finemapped_variant_ld", default = get_or_null(config, "finemapped_variant_ld"))
  parser$add_argument("--disable_sc_cond", default = T, action="store_false", dest="force_sc_replication")

  # Outputs
  parser$add_argument("--output_prefix", default = get_or_null(config, "output_prefix"))

  # Optional overrides
  parser$add_argument("--focal_gene", help = "Gene symbol (e.g. ABCA1)", default = get_or_null(config, "focal_gene"))
  parser$add_argument("--focal_gwas", help = "GWAS trait name or pattern", default = get_or_null(config, "focal_gwas"))

  args <- parser$parse_args(argv)
  output_prefix <- args$output_prefix

  # Now you can use args$hgnc_file, args$focal_gene, etc. as expected
  # ↓↓↓ Call your pipeline using these arguments ↓↓↓

  # ---- Load Data ----
  hgnc <- fread(args$hgnc_file) %>% select(symbol, ensembl_gene_id)
  name_map <- deframe(hgnc %>% select(ensembl_gene_id, symbol))

  cyto <- fread(args$cytoband_file, col.names = c("chr", "start", "end", "band", "v5"))

  cyto_gr <- GRanges(seqnames = cyto$chr,
                     ranges = IRanges(start = cyto$start + 1, end = cyto$end),
                     band = cyto$band)

  variant_annotations <- read_parquet(args$variant_annotations)

  query_gr_variants <- GRanges(seqnames = paste0("chr", variant_annotations$chromosome),
                               ranges = IRanges(start = variant_annotations$bp, end = variant_annotations$bp))
  variant_hits <- findOverlaps(query_gr_variants, cyto_gr)
  variant_annotations$bands <- mcols(cyto_gr[subjectHits(variant_hits)])$band
  variant_annotations$chr_bands <- paste(variant_annotations$chromosome, variant_annotations$bands, sep = "")
  bands_map <- deframe(variant_annotations %>%
                         mutate(loc = paste(chromosome, bands, sep = "")) %>%
                         select(variant, loc))

  gene_reference <- rtracklayer::readGFF(args$gene_reference) %>%
    filter(type == "gene") %>%
    select(gene_id, gene_name, seqid, start, end, strand) %>%
    mutate(TSS = if_else(strand == "+", start, end))

  query_gr_genes <- GRanges(seqnames = paste0("chr", gene_reference$seqid),
                            ranges = IRanges(start = gene_reference$TSS, end = gene_reference$TSS))
  gene_hits <- findOverlaps(query_gr_genes, cyto_gr)
  gene_reference$bands <- NA_character_
  gene_reference$bands[queryHits(gene_hits)] <- mcols(cyto_gr[subjectHits(gene_hits)])$band
  gene_reference$chr_bands <- paste(gene_reference$seqid, gene_reference$bands, sep = "")
  bands_map <- c(bands_map, deframe(gene_reference %>% select(gene_id, chr_bands)))

  cis_trans <- fread(args$cis_trans_file) %>%
    mutate(CisSNPs = str_split(CisSNPs, ";"),
           TransSNPsInStrongLDWithCisSNPs = str_split(TransSNPsInStrongLDWithCisSNPs, ";"))

  gwas_dist_mat <- get_gwas_distance_matrix(args$gwas_dataset, args$uncorrelated_variant_path)

  sc_eqtl_files <- Sys.glob(args$sc_replication)
  force_single_cell_replication <- args$force_sc_replication
  gwas_eqtl_coloc <- fread(args$gwas_eqtl_coloc)

  precalculated_ld_path <- args$finemapped_variant_ld
  ld_matrices <- readRDS(precalculated_ld_path)

  r2_threshold <- 0.5

  variant_index_lookup <- deframe(variant_annotations %>% select(variant, variant_index))

  gwas_eqtl_coloc_filtered <- filter_eqtl_coloc_table(gwas_eqtl_coloc, sc_eqtl_files, force_single_cell_replication)
  gwas_eqtl_coloc_filtered_annot <- gwas_eqtl_coloc_filtered %>%
    mutate(eqtl_variant_index = variant_index_lookup[hit1])

  graph_file_path <- sprintf("%s_graph.rds", output_prefix)
  if (!file.exists(graph_file_path)) {
    g <- build_network(gwas_eqtl_coloc_filtered_annot, cis_trans, variant_annotations, ld_matrices, r2_threshold)
    saveRDS(g, graph_file_path)
  } else {
    g <- readRDS(graph_file_path)
  }

  gene_nodes <- which(V(g)$nodeType == "gene")

  dist_matrix <- get_distance_matrix(g, cached_path = sprintf("%s_distance_matrix.rds", output_prefix))
  path_table <- get_path_table(g, dist_matrix, gene_nodes, cached_path = sprintf("%s_path_table.rds", output_prefix))

  filtered_path_table <- sprintf("%s_path_table_processed.rds", output_prefix)
  if (!file.exists(filtered_path_table)) {
    eqtl_matrix <- as.matrix(fread(args$eqtl_matrix), rownames = "variant")

    path_table_filtered <- filter_path_table(path_table, g, eqtl_matrix) %>%
      as.data.table()

    saveRDS(path_table_filtered, filtered_path_table)
  } else {
    path_table_filtered <- readRDS(filtered_path_table)
  }

  gwas_to_test <- rownames(gwas_dist_mat)

  gene_prioritization_random_variants <- bind_rows(mapply(function(focal_gwas_id, i) {
    message(sprintf("GWAS %s: %s", i, focal_gwas_id))
    random_variant_permutation_test(g, focal_gwas_id, path_table_filtered, gwas_dist_mat)
  }, gwas_to_test, seq_along(gwas_to_test), SIMPLIFY = FALSE), .id = "gwas_id")

  gene_prioritization_random_variants_tibble <- as_tibble(
    gene_prioritization_random_variants %>%
      left_join(hgnc, by = c("gene" = "ensembl_gene_id")) %>%
      arrange(gene_score_p)
  )

  fwrite(gene_prioritization_random_variants_tibble, sprintf("%s_prioritization_perm_results.txt", output_prefix),
         sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  focal_gene <- "BRWD3"
  focal_gwas <- "HEIGHT"

  focal_gene <- args$focal_gene
  focal_gwas <- args$focal_gwas

  if (!is.null(focal_gene) && !is.null(focal_gwas)) {
    focal_gene_id <- hgnc[symbol == focal_gene, ]$ensembl_gene_id
    focal_gwas_id <- grep(focal_gwas, gwas_to_test, value = TRUE)[1]

    if (length(focal_gwas_id) > 0 && length(focal_gene_id) > 0) {
      paths_to_focal_gene <- compute_paths_to_focal_gene(g, path_table_filtered, focal_gwas_id, focal_gene_id)
      plot_focused_network2(g, focal_gwas_id, focal_gene_id, paths_to_focal_gene, name_map, bands_map)
    } else {
      warning("Could not find matching GWAS or gene.")
    }
  }
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}