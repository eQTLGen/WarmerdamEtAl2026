#!/usr/bin/env Rscript
#' This code is responsible for calculating variant level effects with respect to colocalizations to disease
#' and cell type composition traits

# Load libraries
library(tidyverse)
library(data.table)


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

  location <- str_split(variables["available_genes"], "freeze3")[[1]][1]

  available_genes <- fread(variables["available_genes"], header=F)$V1

  gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
    filter(type == "gene") %>% select(gene_id, gene_name, seqid, start, end, strand) %>%
    mutate(TSS = if_else(strand == "+", start, end)) %>%
    filter(gene_id %in% available_genes)

  available_variants <- fread(variables["available_variants"], header=F)$V1

  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type)) %>%
    select(variant, variant_index, phenotype, locus_str, SusieRss_CS, type, typeCisTrans, pip) %>%
    filter(!is.na(pip), phenotype %in% gene_reference$gene_id)

  variant_classes <- fread(variables["variant_coloc_classes"]) %>% select(phenotype, variant_index, variant_type) %>%
    mutate(coloc = variant_type != "Not GWAS colocalising")

  # Process input
  variant_groups <- fread(variables["independent_credible_sets"]) %>%
    separate_wider_delim("gene_id", delim="__", names=c("phenotype", "locus_str", "u_type"), too_many="drop", cols_remove = F) %>%
    mutate(SusieRss_CS = as.numeric(str_extract(cs_id, "\\d+")),
           locus_str = str_replace(locus_str, "_", ":"))

  fine_mapped_variants <- left_join(independent_variant_set_annot_filtered, variant_groups)

  variant_coloc_classes_per_type <- fine_mapped_variants %>%
    inner_join(variant_classes, by = c("variant_index", "phenotype")) %>%
    group_by(type) %>%
    count(coloc) %>%
    mutate(p = n / sum(n) * 100,
           p_lab = sprintf("%s %s (%.1f)", coloc, n, p))

  variant_coloc_classes_per_cs_cluster <- fine_mapped_variants %>%
    inner_join(variant_classes, by = c("variant_index", "phenotype")) %>%
    group_by(CS_group) %>%
    summarise(variant_type_compound = case_when(
      all(variant_type == "Not GWAS colocalising") ~ "Not GWAS colocalising",
      all(variant_type %in% c("Exclusively blood-cell", "Not GWAS colocalising")) ~ "Exclusively blood-cell",
      all(variant_type %in% c("Exclusively non-blood-cell", "Not GWAS colocalising")) ~ "Exclusively non-blood-cell",
      TRUE ~ "Some blood-cell"
    )) %>%
    count(variant_type_compound) %>%
    mutate(p = n / sum(n) * 100,
           p_lab = sprintf("%s %s (%.1f)", variant_type_compound, n, p))

  variant_coloc_classes_per_cs_cluster_per_type <- fine_mapped_variants %>%
    inner_join(variant_classes, by = c("variant_index", "phenotype")) %>%
    group_by(type, CS_group) %>%
    summarise(variant_type_compound = case_when(
      all(variant_type == "Not GWAS colocalising") ~ "Not GWAS colocalising",
      all(variant_type %in% c("Exclusively blood-cell", "Not GWAS colocalising")) ~ "Exclusively blood-cell",
      all(variant_type %in% c("Exclusively non-blood-cell", "Not GWAS colocalising")) ~ "Exclusively non-blood-cell",
      TRUE ~ "Some blood-cell"
    )) %>%
    count(variant_type_compound) %>%
    mutate(p = n / sum(n) * 100,
           p_lab = sprintf("%s %s (%.1f)", variant_type_compound, n, p))

  ggplot(variant_coloc_classes_per_cs_cluster_per_type, aes(x="", y=p, fill=variant_type_compound)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    geom_text_repel(aes(x =  1.2,label = p_lab), position = position_stack(vjust = 0.5)) +
    facet_wrap(vars(type), scales="free")

  # ---
  # Count number of eQTLs for which we identify GWAS colocalizations
  coloc_glob <- variables["coloc_cleaned"]

  colocalizations <- bind_rows(mapply(fread, Sys.glob(coloc_glob), SIMPLIFY = F)) %>%
    filter(eQTL_gene %in% gene_reference$gene_id) %>%
    mutate(gwas_coloc = T)

  independent_variant_set_annot_filtered_flagged <- left_join(
    independent_variant_set_annot_filtered,
    colocalizations %>% distinct(eQTL_gene, signal1, eQTL_locus, gwas_coloc),
    by = c("phenotype" = "eQTL_gene", "SusieRss_CS" = "signal1", "locus_str" = "eQTL_locus")) %>%
    left_join(variant_groups)

  independent_variant_set_annot_filtered_flagged %>%
    summarise(n_all = n(), n_gwas_coloc = sum(gwas_coloc, na.rm=T), p = n_gwas_coloc / n_all * 100)

  independent_variant_set_annot_filtered_flagged %>%
    group_by(typeCisTrans) %>%
    summarise(n_all = n(), n_gwas_coloc = sum(gwas_coloc, na.rm=T), p = n_gwas_coloc / n_all * 100)

  independent_variant_set_annot_filtered_flagged %>%
    group_by(typeCisTrans, CS_group) %>%
    summarise(gwas_coloc = any(gwas_coloc, na.rm=T)) %>%
    group_by(typeCisTrans) %>%
    summarise(n_all = n(), n_gwas_coloc = sum(gwas_coloc, na.rm=T), p = n_gwas_coloc / n_all * 100)
  # Work
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}