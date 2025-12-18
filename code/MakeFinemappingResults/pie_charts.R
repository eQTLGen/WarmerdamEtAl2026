#!/usr/bin/env Rscript


# Load libraries
library(data.table)
library(arrow)
library(tidyverse)
library(dtplyr)
library(dbplyr)
library(igraph)
library(ggrepel)

source("code/ComparisonPhase1/compare_against_phase1.R")

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

  ordered_variant_types <- c("Exclusively blood-cell",
                             "Some blood-cell",
                             "Exclusively non-blood-cell",
                             "Not GWAS colocalising")

  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))
  location <- str_split(variables["available_genes"], "freeze3")[[1]][1]

  available_genes <- fread(variables["available_genes"], header=F)$V1
  outcome_variables["n_available_genes_meta_analysis"] <- length(available_genes)

  available_variants <- fread(variables["available_variants"], header=F)$V1

  variant_classes <- fread(variables["variant_coloc_classes"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type))

    class_counts_eqtl_level <- variant_classes %>%
    filter(!is.na(pip), !is.na(variant_type)) %>%
    mutate(variant_type = factor(variant_type, levels = ordered_variant_types, ordered = T)) %>%
    group_by(type) %>%
    count(variant_type) %>%
    mutate(p = n / sum(n) * 100,
           p_lab = sprintf("%s %s (%.1f)", variant_type, n, p))

  coloc_counts_eqtl_level <- variant_classes %>%
    filter(!is.na(pip), !is.na(variant_type)) %>%
    mutate(coloc = variant_type != "Not GWAS colocalising") %>%
    group_by(type != "intermediate") %>%
    count(coloc) %>%
    mutate(p = n / sum(n) * 100,
           p_lab = sprintf("%s %s (%.1f)", coloc, n, p))


  ggplot(class_counts_eqtl_level, aes(x="", y=p, fill=variant_type)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    geom_text_repel(aes(x =  1.2,label = p_lab), position = position_stack(vjust = 0.5)) +
    facet_wrap(vars(type), scales="free")

  cis_trans_relationships <- fread(variables["cis_trans_relationships"])

  # For every trans-eQTL variant, identify if it colocalizes to a cis-eGene, and for how many CSs this is the case
  cis_trans_coloc <- fread(variables["lbf_based_cis_trans_coloc"])

  variant_annotations <- fread("independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-maxR2_0.9-noHla-noCrossmapping-annotations-v12.txt.gz")

  trans_gene_pairs_with_cis_coloc <- cis_trans_coloc %>%
    inner_join(cis_trans_relationships, by = c("cis_gene" = "CisGene", "trans_gene" = "TransGene")) %>%
    group_by(trans_gene, variant_index_hit2) %>%
    summarise(MaxNrCredibleSetsAffectingBothCisAndTrans = max(NrCredibleSetsAffectingBothCisAndTrans, na.rm=T))

  View(trans_gene_pairs_with_cis_coloc)

  class_counts_eqtl_level_stratified_by_nr_of_colocalizing_cis_signals <- variant_classes_with_cis_trans_relationships %>%
    filter(type %in% c("trans"), !is.na(pip)) %>%
    mutate(variant_type = factor(variant_type, levels = ordered_variant_types, ordered = T)) %>%
    group_by(type, MaxNrCredibleSetsAffectingBothCisAndTrans) %>%
    count(variant_type) %>%
    mutate(p = n / sum(n) * 100,
           p_lab = sprintf("%s %s (%.1f)", variant_type, n, p))

  ggplot(class_counts_eqtl_level_stratified_by_nr_of_colocalizing_cis_signals, aes(x="", y=p, fill=variant_type)) +
    geom_bar(stat="identity", width=1) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
    facet_grid(rows=vars(type), cols=vars(MaxNrCredibleSetsAffectingBothCisAndTrans), scales="free")

  ggsave("trans-eQTLs_stratified_across_cis-eQTL_colocalization.pdf")


  variant_annotations_with_cis_trans_relationships <- left_join(
    variant_annotations, trans_gene_pairs_with_cis_coloc,
    by=c("phenotype" = "trans_gene", "variant_index" = "variant_index_hit2")) %>%
    mutate(MaxNrCredibleSetsAffectingBothCisAndTrans = replace_na(MaxNrCredibleSetsAffectingBothCisAndTrans, 0))


  rbs_eqtl_level_stratified_by_nr_of_colocalizing_cis_signals <- variant_annotations_with_cis_trans_relationships %>%
    filter(type %in% c("trans"), !is.na(pip)) %>%
    group_by(type, MaxNrCredibleSetsAffectingBothCisAndTrans) %>%
    filter(!is.na(scMeta_CD4_T_beta), !is.na(scMeta_CD4_T_se)) %>%
    summarise(rb_mat = calcu_cor_true(beta, standard_error, scMeta_CD4_T_beta, scMeta_CD4_T_se, theta = 1)) %>%
    unnest(cols = c(rb_mat), names_sep="_") %>%
    mutate(rb_lab = sprintf("Rb of betas: %.3f (SE = %.3f)", rb_mat_r, rb_mat_se_r))

  ggplot(rbs_eqtl_level_stratified_by_nr_of_colocalizing_cis_signals, aes(x=rb_mat_r, y=MaxNrCredibleSetsAffectingBothCisAndTrans)) +
    geom_point() +
    geom_errorbarh(aes(xmin = rb_mat_r-rb_mat_se_r, xmax = rb_mat_r+rb_mat_se_r), height = 0.2, color = "gray40") +
    geom_vline(xintercept=0)


    ggsave("trans-eQTL_rbs_inCd4TCells_across_cis-eQTL_colocalization.pdf")

  variant_classes_of_cis_eqtls_with_trans_relationships <- left_join(
    variant_classes, cis_trans_coloc,
    by=c("phenotype" = "cis_gene", "variant_index" = "variant_index"))

  class_counts_eqtl_level_stratified_by_colocalizing_cis_signals <- variant_classes_of_cis_eqtls_with_trans_relationships %>%
    filter(type %in% c("cis"), !is.na(pip)) %>%
    mutate(variant_type = factor(variant_type, levels = ordered_variant_types, ordered = T),
           trans_effect = !is.na(trans_gene)) %>%
    group_by(type, trans_effect) %>%
    count(variant_type) %>%
    mutate(p = n / sum(n) * 100,
           p_lab = sprintf("%s %s (%.1f)", variant_type, n, p))

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}