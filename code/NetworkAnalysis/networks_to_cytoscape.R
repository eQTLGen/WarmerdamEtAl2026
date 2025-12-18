#!/usr/bin/env Rscript


# Load libraries
library(igraph)
library(data.table)
library(arrow)
library(tidyverse)
library(GenomicRanges)
library(ggraph)
library(RCy3)
library(RColorBrewer)
library(broom)
library(tidyr)

# Declare constants
non_blood_cell_type_traits <- c("PASS_CD_deLange2017", "PASS_Celiac", "PASS_IBD_deLange2017",
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


ideal_rect_size <- function(text,
                            family="HelveticaNeue", size = 12,
                            padding_mm = 2) {

  # Measure text dimensions
  w_in <- strwidth(text, family = family, units = "inch", font=size)
  h_in <- strheight(text, family = family, units = "inch", font=size)

  # Padding
  padding_in <- padding_mm / 25.4
  w_in <- w_in + 2 * padding_in
  h_in <- h_in + 2 * padding_in

  # Convert inches -> Cytoscape units (96 dpi = 96 pixels/inch)
  w_px <- w_in * 72
  h_px <- h_in * 72

  return(list(width_pt = w_px, height_pt = h_px))
}

# Declare function definitions
set_style_lh <- function() {
  # --------------------------------------------------------
  # Define Cytoscape style using mapping functions
  # --------------------------------------------------------
  style.name <- "grn_style_lh"
  if (style.name %in% getVisualStyleNames()) {
    deleteVisualStyle(style.name)
  }
  copyVisualStyle("default", style.name)

  # Node mappings
  setNodeShapeDefault("ELLIPSE", style.name)
  setNodeLabelMapping("label",
                      style.name = style.name)

  setNodeBorderColorMapping('node_type',
                            c("focal","gene","association"),
                            c("black","black","#666666"),
                            mapping.type = 'd',
                            style.name = style.name)

  setNodeBorderWidthMapping('node_type',
                            c("focal","gene","association"),
                            c(2, 2, 0),
                            mapping.type = 'd',
                            style.name = style.name)

  setNodeColorMapping('node_type',
                          c("focal","gene","association"),
                          c("white","white","#B3B3B3"),
                      mapping.type = 'd',
                      style.name = style.name)

  setNodeSizeMapping('node_type',
                     c("focal","gene","association"),
                     c(50,30,20),
                     mapping.type = 'd',
                     style.name = style.name)


  # Edge mappings
  setEdgeColorMapping('Direction',
                      c("1","-1"),
                      c("#5AB847","#E41E25"),
                      mapping.type = 'd',
                      style.name = style.name)

  setEdgeLineStyleMapping('edgeType',
                          c("genepairs","cis","trans"),
                          c("SOLID","SOLID","LONG_DASH"),
                          style.name = style.name)

  setEdgeTargetArrowShapeMapping('Direction',
                                 c("1","-1"),
                                 c("DELTA_SHORT_1","T"),
                                 style.name = style.name)

  setEdgeTargetArrowColorMapping('Direction',
                                 c("1","-1"),
                                 c("#5AB847","#E41E25"),
                                 mapping.type = 'd',
                                 style.name = style.name)

  setVisualStyle(style.name)
}


# Declare function definitions
define_style_p2 <- function() {
  # --------------------------------------------------------
  # Define Cytoscape style using mapping functions
  # --------------------------------------------------------
  style.name <- "grn_style_p2"
  if (style.name %in% getVisualStyleNames()) {
    deleteVisualStyle(style.name)
  }
  copyVisualStyle("default", style.name)


  setNodeFontFaceDefault(c("HelveticaNeue"), style.name=style.name)
  setEdgeFontFaceDefault(c("HelveticaNeue"), style.name=style.name)

  # Node mappings
  setNodeLabelMapping("label",
                      style.name = style.name)

  setNodeFontFaceMapping('node_type',
                         c("focal","gene","association"),
                         c("HelveticaNeue-MediumItalic", "HelveticaNeue-Italic", "HelveticaNeue"),
                         style.name=style.name)

  setNodeBorderColorMapping('node_type',
                            c("focal","gene","association"),
                            c("#000000","#000000","#666666"),
                            mapping.type = 'd',
                            style.name = style.name)

  setNodeBorderWidthMapping('node_type',
                            c("focal","gene","association"),
                            c(3, 3, 3),
                            mapping.type = 'd',
                            style.name = style.name)

  setNodeShapeMapping('node_type',
                           c("focal","gene","association"),
                           c("ROUND_RECTANGLE", "ROUND_RECTANGLE", "ELLIPSE"),
                           style.name = style.name)

  setNodeLabelColorMapping('node_type',
                           c("focal","gene","association"),
                           c("black", "black", "black"),
                           mapping.type = 'd',
                           style.name = style.name)

  setNodeColorMapping('node_type',
                      c("focal","gene","association"),
                      c("#ffffff","#ffffff","#B3B3B3"),
                      mapping.type = 'd',
                      style.name = style.name)

  setNodeWidthMapping('node_width',
                     mapping.type = 'p',
                     style.name = style.name)

  setNodeHeightMapping('node_height',
                      mapping.type = 'p',
                      style.name = style.name)


  # Edge mappings
  setEdgeColorMapping('Direction',
                      c("1","-1"),
                      c("#5AB847","#E41E25"),
                      mapping.type = 'd',
                      style.name = style.name)

  setEdgeLineStyleMapping('edge_type',
                          c("genepairs","cis","trans"),
                          c("SOLID","SOLID","EQUAL_DASH"),
                          style.name = style.name)

  setEdgeTargetArrowShapeMapping('Direction',
                                 c("1","-1"),
                                 c("DELTA_SHORT_1","T"),
                                 style.name = style.name)

  setEdgeTargetArrowColorMapping('Direction',
                                 c("1","-1"),
                                 c("#5AB847","#E41E25"),
                                 mapping.type = 'd',
                                 style.name = style.name)

  setEdgeLabelMapping('sc_replication',
                      style.name = style.name)
}

set_style_p2 <- function() {
  style.name <- "grn_style_p2"
  setVisualStyle(style.name)
}


get_focused_network <- function(g, focal_gene_id,
                                     paths_to_focal_gene, hgnc_name_map,
                                     bands_map) {

  # --------------------------------------------------------
  # Subgraph with relevant nodes
  # --------------------------------------------------------
  relevant_nodes <- na.omit(c(paths_to_focal_gene$source,
                              paths_to_focal_gene$intermediate_node,
                              focal_gene_id))
  g_sub <- induced_subgraph(g, vids = relevant_nodes)

  # --------------------------------------------------------
  # Build node table
  # --------------------------------------------------------

  nodes <- igraph::as_data_frame(g_sub, "vertices") %>%
    mutate(
      id = name,
      # Assign node types
      node_type = case_when(
        name == focal_gene_id ~ "focal",
        nodeType == "gene" ~ "gene",
        nodeType == "association" ~ "association"
      ),

      # Assign bands
      bands = case_when(
        nodeType == "gene" ~ bands_map[name],
        nodeType == "association" ~ bands_map[index_variant]
      ),

      # Assign node labels
      label = unname(case_when(
        name %in% names(hgnc_name_map) ~ sprintf("%s", hgnc_name_map[name]),
        node_type %in% c("association") & !is.na(index_variant) ~ bands_map[index_variant],
        node_type == "association" ~ bands_map[str_extract(name, "(\\w+):(\\w+)", group=2)]
      ))
    )

  rect_sizes <- ideal_rect_size(nodes$label)

  nodes$node_width <- if_else(nodes$node_type == "association", 8, rect_sizes$width_pt)
  nodes$node_height <- if_else(nodes$node_type == "association", 8, rect_sizes$height_pt)

  edges <- igraph::as_data_frame(g_sub, "edges") %>%
    mutate(
      source = from,
      target = to,
      # Assign edge types and properties
      edge_type = case_when(
        edgeType == "genelink" ~ "genepairs",
        typeCisTrans == "cis" ~ "cis",
        typeCisTrans == "trans" ~ "trans"
      ),
      edge_sign = Direction,
      edge_weight = if_else(edge_type == "ld", 1, 10),
      arrow_type = "closed",
      sc_sum_of_z = rowSums(pick(starts_with("sc_z_score"))),
      sc_z_score_meta = sc_sum_of_z / sqrt(ncol(pick(starts_with("sc_z_score")))),
      sc_replication = if_else(is.na(cell_type), "", sprintf("%.2f (%s)", risk_allele_sc_z, cell_type)),
      sc_replication_meta = if_else(is.na(cell_type), "", sprintf("%.2f", sc_z_score_meta))
    )

  return(list(nodes = nodes, edges = edges))
}


# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  source("code/NetworkAnalysis/network_analyses_coloc.R")
  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))

  prefix <- variables["network_analyses_prefix"]

  hgnc <- fread(variables["hgnc_complete_set"]) %>% select(ensembl_gene_id, symbol)
  hgnc_name_map <- deframe(hgnc)

  variant_annotations <- read_parquet(variables["variant_annotations"])

  cyto <- fread(variables["cytoband_reference"], col.names = c("chr", "start", "end", "band", "v5"))

  cyto_gr <- GRanges(seqnames = cyto$chr,
                     ranges = IRanges(start = cyto$start + 1, end = cyto$end),
                     band = cyto$band)

  query_gr_variants <- GRanges(seqnames = paste0("chr", variant_annotations$chromosome),
                               ranges = IRanges(start = variant_annotations$bp, end = variant_annotations$bp))
  variant_hits <- findOverlaps(query_gr_variants, cyto_gr)
  variant_annotations$bands <- mcols(cyto_gr[subjectHits(variant_hits)])$band
  variant_annotations$chr_bands <- paste(variant_annotations$chromosome, variant_annotations$bands, sep = "")

  bands_map <- deframe(variant_annotations %>%
                         mutate(loc = paste(chromosome, bands, sep = "")) %>%
                         select(variant, loc))

  g <- readRDS(paste0(prefix, "_graph.rds"))

  association_nodes_df <- tibble(igraph::as_data_frame(g, what = "vertices")[V(g)$nodeType == "association", ]) %>%
    select(name, gwas_id)

  path_table_filtered <- readRDS(paste0(prefix, "_path_table_processed.rds")) %>%
    inner_join(association_nodes_df, by = c("association" = "name"))

  cytoscapePing ()
  cytoscapeVersionInfo ()

  define_style_p2()

  focal_gene <- "ABCA1"
  focal_gwas <- "Low"

  focal_gene_id <- hgnc[symbol == focal_gene, ]$ensembl_gene_id
  focal_gwas_id <- grep(focal_gwas, non_blood_cell_type_traits, value = TRUE, ignore.case = TRUE)[1]

  print(focal_gene_id)
  print(focal_gwas_id)

  if (length(focal_gwas_id) > 0 && length(focal_gene_id) > 0) {
    paths_to_focal_gene <- compute_paths_to_focal_gene(g, path_table_filtered, focal_gwas_id, focal_gene_id)

    node_edge_list <- get_focused_network(g=g, focal_gene_id = focal_gene_id,
                             paths_to_focal_gene = paths_to_focal_gene, hgnc = hgnc_name_map, bands_map = bands_map)

    title_text <- sprintf("%s - %s", focal_gwas_id, focal_gene_id)

    createNetworkFromDataFrames(nodes = node_edge_list$nodes, edges = node_edge_list$edges, title= title_text,
                                collection = "Focused_Networks")

    set_style_p2()
  } else {
    warning("Could not find matching GWAS or gene.")
  }
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}