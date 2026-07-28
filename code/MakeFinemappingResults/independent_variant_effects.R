#!/usr/bin/env Rscript
#' This code is responsible for calculating variant level effects with respect to colocalizations to disease
#' and cell type composition traits

# Load libraries
library(tidyverse)
library(data.table)
library(ggrepel)
library(ggideogram)
library(ggpubr)


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

  hgnc_symbols <- fread(variables["hgnc_symbols"]) %>% select(ensembl_gene_id, symbol) %>% deframe()

  available_genes <- fread(variables["available_genes"], header=F)$V1

  gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
    filter(type == "gene") %>% select(gene_id, gene_name, seqid, start, end, strand) %>%
    mutate(TSS = if_else(strand == "+", start, end)) %>%
    filter(gene_id %in% available_genes)

  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"]) %>%
    mutate(typeCisTrans = if_else(type == "intermediate", "trans", type)) %>%
    select(variant, variant_index, phenotype, locus_str, SusieRss_CS, type, typeCisTrans, pip, chromosome, bp) %>%
    filter(phenotype %in% gene_reference$gene_id)

  variant_classes <- fread(variables["variant_coloc_classes"]) %>% select(phenotype, variant_index, variant_type) %>%
    mutate(coloc = variant_type != "Not GWAS colocalising")

  # Process input
  variant_groups <- fread(variables["independent_credible_sets"]) %>%
    separate_wider_delim("gene_id", delim="__", names=c("phenotype", "locus_str", "u_type"), too_many="drop", cols_remove = F) %>%
    mutate(SusieRss_CS = as.numeric(str_extract(cs_id, "\\d+")),
           locus_str = str_replace(locus_str, "_", ":"))

  fine_mapped_variants <- left_join(independent_variant_set_annot_filtered, variant_groups) %>%
    filter(!is.na(SusieRss_CS))

  variant_coloc_classes_per_type <- fine_mapped_variants %>%
    inner_join(variant_classes, by = c("variant_index", "phenotype")) %>%
    group_by(type) %>%
    count(coloc) %>%
    mutate(p = n / sum(n) * 100,
           p_lab = sprintf("%s %s (%.1f)", coloc, n, p))

  cis_trans_coloc_proc <- fread(variables["lbf_based_cis_trans_coloc"]) %>%
    select(cis_gene = "cis_gene", trans_gene = "trans_gene", cis_variant = "hit1", trans_variant = "hit2") %>%
    mutate(cis_symbol = hgnc_symbols[cis_gene])

  trans_pair_annotations <- cis_trans_coloc_proc %>%
    group_by(trans_variant, trans_gene) %>%
    summarise(cis_symbols = list(unique(cis_symbol)))

  variant_coloc_classes_per_cs_cluster <- fine_mapped_variants %>%
    inner_join(variant_classes, by = c("variant_index", "phenotype")) %>%
    group_by(CS_group, type == "cis") %>%
    summarise(variant_type_compound = case_when(
      all(variant_type == "Not GWAS colocalising") ~ "Not GWAS colocalising",
      all(variant_type %in% c("Exclusively blood-cell", "Not GWAS colocalising")) ~ "Exclusively blood-cell",
      all(variant_type %in% c("Exclusively non-blood-cell", "Not GWAS colocalising")) ~ "Exclusively non-blood-cell",
      TRUE ~ "Some blood-cell"
    )) %>%
    count(variant_type_compound) %>%
    mutate(p = n / sum(n) * 100,
           p_lab = sprintf("%s %s (%.1f)", variant_type_compound, n, p))

  cis_gene_limit <- 5

  get_cis_gene_lab <- function(cis_genes, cis_genes_counts, n_trans_genes) {

    n_cis_genes <- length(cis_genes)
    n_cis_genes_to_show <- ifelse(n_cis_genes <= cis_gene_limit, n_cis_genes, cis_gene_limit - 1)

    if(is.na(n_cis_genes) | is.null(n_cis_genes) | n_cis_genes == 0) {
      return("None")
    }

    gene_list <- paste0("", cis_genes[1:n_cis_genes_to_show], " (", round(cis_genes_counts[1:n_cis_genes_to_show] / n_trans_genes * 100), "%)", collapse=",\n")

    extra_genes <- if_else(
      n_cis_genes <= cis_gene_limit, "",
      sprintf(",\n+%s more", n_cis_genes - n_cis_genes_to_show))

    full_lab <- sprintf(
      fmt= "%s%s", gene_list, extra_genes)

    return(full_lab)
  }

  variant_coloc_classes_per_cs_cluster_per_type <- fine_mapped_variants %>%
    inner_join(variant_classes, by = c("variant_index", "phenotype")) %>%
    left_join(trans_pair_annotations, by = c("variant" = "trans_variant", "phenotype" = "trans_gene")) %>%
    group_by(chromosome, type, CS_group) %>%
    summarise(
      variant_type_compound = case_when(
        all(variant_type == "Not GWAS colocalising") ~ "Not GWAS colocalising",
        all(variant_type %in% c("Exclusively blood-cell", "Not GWAS colocalising")) ~ "Exclusively blood-cell",
        all(variant_type %in% c("Exclusively non-blood-cell", "Not GWAS colocalising")) ~ "Exclusively non-blood-cell",
        TRUE ~ "Some blood-cell"
      ),
      n_genes = n(),
      cis_genes_counts = list(sort(table(unlist(cis_symbols)), decreasing = TRUE)),
      cis_genes = list(names(sort(table(unlist(cis_symbols)), decreasing=T))),
      cis_genes_lab = get_cis_gene_lab(unlist(cis_genes), unlist(cis_genes_counts), n_genes),
      bp_min = min(bp),
      bp_max = max(bp),
      bp_mean = mean(bp),
    ) %>% ungroup()

  variant_coloc_classes_per_cs_cluster_per_type_summarised <- variant_coloc_classes_per_cs_cluster_per_type %>%
    mutate(typeCisTrans = if_else(type == "cis", "cis", "trans")) %>%
    group_by(variant_type_compound, type) %>% summarise(n = n(), n_genes = sum(n_genes)) %>%
    group_by(type) %>%
    mutate(p = n / sum(n) * 100,
           p_genes = n_genes / sum(n_genes) * 100,
           p_lab = sprintf("%s %s (%.1f)", variant_type_compound, n, p))

  variant_coloc_classes_per_cs_cluster_per_type %>%
    filter(type == "trans") %>%
    ungroup() %>%
    summarise(
      n_equals_1 = sum(n_genes == 1),
      n_total = n(),
      n_over_100 = sum(n_genes > 100),
      n_equals_1_p = n_equals_1 / n_total * 100,
      n_over_100_p = n_over_100 / n_total * 100)

  trans_only_variant_annots <- variant_coloc_classes_per_cs_cluster_per_type %>% filter(type == "trans") %>%
    mutate(trans_only = if_else(cis_genes_lab == "None", "trans-only", "pair"),
           gwas_coloc = if_else(variant_type_compound != "Not GWAS colocalising", "GWAS colocalizing", "Not GWAS colocalizing"),
           bctc_coloc = variant_type_compound %in% c("Exclusively blood-cell", "Some blood-cell")
    )

  freq_table <- table(trans_only_variant_annots$trans_only,
                      trans_only_variant_annots$gwas_coloc)
  freq_table
  fisher_res <- fisher.test(freq_table)
  fisher_res
  fisher_res$p.value

  freq_table_no_bctc <- table(trans_only_variant_annots$trans_only[!trans_only_variant_annots$bctc_coloc],
                              trans_only_variant_annots$gwas_coloc[!trans_only_variant_annots$bctc_coloc])
  freq_table_no_bctc
  fisher_res_no_bctc <- fisher.test(freq_table_no_bctc)
  fisher_res_no_bctc
  fisher_res_no_bctc$p.value

  variant_coloc_classes_per_cs_cluster_per_type_summarised <- variant_coloc_classes_per_cs_cluster_per_type %>%
    count(variant_type_compound, type) %>%
    group_by(type) %>%
    mutate(p = n / sum(n) * 100,
           p_lab = sprintf("%s %s (%.1f)", variant_type_compound, n, p))

  ggplot(variant_coloc_classes_per_cs_cluster_per_type_summarised, aes(x="", y=p, fill=variant_type_compound)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    geom_text_repel(aes(x =  1.2,label = p_lab), position = position_stack(vjust = 0.5)) +
    facet_wrap(vars(type), scales="free")

  ggsave("variant_level_pie_charts.pdf", width=11, height=11)

  # ---
  # Count number of eQTLs for which we identify GWAS colocalizations
  coloc_glob <- variables["coloc_cleaned"]

  colocalizations <- bind_rows(mapply(fread, Sys.glob(coloc_glob), SIMPLIFY = F)) %>%
    filter(eQTL_gene %in% gene_reference$gene_id) %>%
    mutate(gwas_coloc = T)

  independent_variant_set_annot_filtered_flagged <- left_join(
    independent_variant_set_annot_filtered %>% filter(!is.na(SusieRss_CS)),
    colocalizations %>% distinct(eQTL_gene, signal1, eQTL_locus, gwas_coloc),
    by = c("phenotype" = "eQTL_gene", "SusieRss_CS" = "signal1", "locus_str" = "eQTL_locus")) %>%
    left_join(variant_groups)

  independent_variant_set_annot_filtered_flagged %>%
    summarise(n_all = n(), n_gwas_coloc = sum(gwas_coloc, na.rm=T), p = n_gwas_coloc / n_all * 100)

  independent_variant_set_annot_filtered_flagged %>%
    group_by(typeCisTrans) %>%
    summarise(n_all = n(), n_gwas_coloc = sum(gwas_coloc, na.rm=T), p = n_gwas_coloc / n_all * 100)

  independent_variant_set_annot_filtered_flagged %>%
    group_by(type, CS_group) %>%
    summarise(gwas_coloc = any(gwas_coloc, na.rm=T)) %>%
    group_by(type) %>%
    summarise(n_all = n(), n_gwas_coloc = sum(gwas_coloc, na.rm=T), p = n_gwas_coloc / n_all * 100)
  # Work
  # Process output


  RollingSum <- function(x, window = 1e+06) {
    Sum <- rep(0, times = length(x))
    for (i in 1:length(x)) {
      Sum[i] <- length(x[x > x[i] - window & x < x[i] +
        window])
    }
    return(Sum)
  }

  window <- 1e6

  independent_variant_set_annot_filtered_rolling_sum <- independent_variant_set_annot_filtered %>%
    group_by(chromosome) %>%
    arrange(bp) %>%
    mutate(
      count_1Mb = sapply(bp, function(x) {
        sum(bp >= x - window/2 & bp <= x + window/2)
      })
    ) %>%
    ungroup()

  coloc_colors <- c("Exclusively non-blood-cell" = "#44a7dd",
                    "Exclusively blood-cell" = "#b66e11",
                    "Some blood-cell" = "#fcd92e",
                    "Not GWAS colocalising" = "grey70")

  chromosome_row_col_def <- tibble(chromosome = 1:22) %>%
    mutate(chrom_row = (chromosome-1) %/% 5, chrom_col = chromosome - (chrom_row * 5))

  data(karyotype_hg38, package="ggideogram")

  karyotype_hg38_row_col <- karyotype_hg38 %>%
    mutate(chromosome = as.integer(str_remove(Chrom, "chr"))) %>%
    inner_join(chromosome_row_col_def)

  hotspots <- ggplot(data = variant_coloc_classes_per_cs_cluster_per_type %>%
    filter(type != "cis") %>% inner_join(chromosome_row_col_def)) +
    geom_ideogram(data = karyotype_hg38_row_col, aes(x = -100, ymin = End, ymax = Start,
                                                     chrom = Chrom, fill = Stain, arm = Arm),
                  radius = unit(1, 'pt'), width = 80,
                  linewidth = (0.25 / (ggplot2::.pt * 72.27/96)), colour = 'black',
                  show.legend = FALSE) +
    scale_fill_manual(values = cytoband_colors) +
    geom_segment(aes(y = bp_mean, yend = bp_mean, x=0, xend=n_genes, color = variant_type_compound),
                 inherit.aes = F, linewidth = (1 / (ggplot2::.pt * 72.27/96))) +
    scale_y_reverse(name = "Location on chromosome in million basepairs (hg38)", labels = ~ sprintf("%.1f",.x/1e6)) +
    scale_x_continuous(name = "Number of genes associated")+
    coord_cartesian(xlim = c(-150, 1500), clip = "off") +
    scale_color_manual(values = coloc_colors) +
    geom_text_repel(data = . %>% filter(n_genes > 100), aes(y=bp_mean, x=n_genes, label=cis_genes_lab),
                    size=4.5/.pt, hjust = 0, nudge_x = 100, nudge_y=-5e6, max.overlaps = Inf, min.segment.length = 0, parse=F) +
    facet_grid(rows = vars(chrom_row), cols = vars(chrom_col), scales = "free_y", space="free_y") +
    theme(legend.position = "bottom")

  ggsave("genome_hotspots_trans.pdf", plot=hotspots, width=180, height=240, units="mm")

  df_binned <- variant_coloc_classes_per_cs_cluster_per_type %>%
    filter(type == "trans") %>%
    arrange(n_genes) %>%
    mutate(bin = case_when(
      n_genes == 1 ~ "1",
      n_genes == 2 ~ "2",
      n_genes == 3 ~ "3",
      n_genes == 4 ~ "4",
      n_genes >= 5 & n_genes <= 9 ~ "5-9",
      n_genes >= 10 & n_genes <= 19 ~ "10-19",
      n_genes >= 20 & n_genes <= 49 ~ "20-49",
      n_genes >= 50 & n_genes <= 99 ~ "50-99",
      n_genes >= 100 & n_genes <= 199 ~ "100-199",
      n_genes >= 200 & n_genes <= 499 ~ "200-499",
      n_genes >= 500 & n_genes <= 999 ~ "500-999",
      n_genes >= 1000 ~ "1000+",
      TRUE ~ NA_character_
    ), bin = ordered(bin, levels=unique(bin)))

  hotspots2 <- ggplot(data = df_binned, aes(bin)) + geom_bar() +
    xlab("Number of trans-eQTLs per independent genetic signal") +
    ylab("Count of independent genetic signals") + ggtitle("Trans-eQTL hotspots")

  ggsave("trans_eQTL_hotspots_hist.pdf", plot=hotspots2, width=120, height=90, units="mm")
  ggsave("trans_eQTL_hotspots_hist.png", plot=hotspots2, width=120, height=90, units="mm")

  df_binned_summarised <- df_binned %>%
    group_by(bin, variant_type_compound) %>%
    tally() %>% group_by(variant_type_compound) %>% mutate(p = n / sum(n))
  hotspots3 <- ggplot(df_binned_summarised, aes(bin, y=p, fill=variant_type_compound)) +
    geom_col(position = position_dodge()) +
    scale_color_manual(values = coloc_colors) +
    xlab("Number of trans-eQTLs per independent genetic signal") +
    ylab("Count of independent genetic signals") + ggtitle("Trans-eQTL hotspots")

  hotspots3

  df_binned_5_genes <- df_binned %>% mutate(over_five_trans_eqtls = n_genes > 5,
                                            blood_cell = variant_type_compound == "Exclusively blood-cell") %>%
    mutate(variant_type_compound_fac = ordered(variant_type_compound, levels=c("Not GWAS colocalising", "Exclusively non-blood-cell", "Some blood-cell", "Exclusively blood-cell")))

  df_binned_summarised_2 <- df_binned_5_genes %>%
    group_by(over_five_trans_eqtls, blood_cell) %>%
    tally() %>% group_by(blood_cell) %>% mutate(p = n / sum(n))

  df_binned_5_genes_summarised <- df_binned_5_genes %>%
    group_by(variant_type_compound_fac) %>%
    summarise(summ = list(summary(n_genes))) %>% unnest_wider(summ)

  comparisons <- list(
    c("Exclusively blood-cell", "Exclusively non-blood-cell"),
    c("Exclusively blood-cell", "Not GWAS colocalising"),
    c("Exclusively blood-cell", "Some blood-cell"),
    c("Exclusively non-blood-cell", "Some blood-cell"),
    c("Some blood-cell", "Not GWAS colocalising"),
    c("Not GWAS colocalising", "Exclusively non-blood-cell")
  )

  library(rstatix)

  stat.test <- df_binned_5_genes %>%
    pairwise_wilcox_test(
      n_genes ~ variant_type_compound_fac,
      comparisons = comparisons
    )
  stat.test$label <- sprintf("%.2e", stat.test$p)
  stat.test$y.position <- c(3.8, 3.6, 4.0, 3.2, 3.4, 3.0)

  coloc_colors_bright <- c("Exclusively non-blood-cell" = "#44a7dd",
                           "Exclusively blood-cell" = "#f5911e",
                           "Some blood-cell" = "#fcd92e",
                           "Not GWAS colocalising" = "#c9cbcd")

  coloc_outline_bright <- c("Exclusively non-blood-cell" = "#0056a0",
                            "Exclusively blood-cell" = "#b66e11",
                            "Some blood-cell" = "#d5ab2c",
                            "Not GWAS colocalising" = "grey70")


  trans_eqlts_between_classes <- ggplot(df_binned_5_genes, aes(variant_type_compound_fac, n_genes, fill =variant_type_compound_fac, color =variant_type_compound_fac)) +
    geom_violin() +
    scale_fill_manual(values=coloc_colors_bright) +
    scale_color_manual(values=coloc_outline_bright) +
    geom_boxplot(width = 0.05, fill="white") +
    scale_y_log10(name="Number of emanating trans-eQTLs") +
    geom_text(data=df_binned_5_genes_summarised,
              aes(y = 0.6, label=sprintf("Median = %.0f,\nMean = %.2f", Median, Mean)),
              size=6/.pt) +
    stat_pvalue_manual(
      stat.test,
      label = "p",
      size=6/.pt,
      p.format.style = "scientific"
    ) + theme(legend.position = "none")

  trans_eqlts_between_classes

  ggsave("trans_eQTLs_between_classes.pdf", width=150, height=120, units="mm")

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}