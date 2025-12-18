#!/usr/bin/env Rscript
#' This code recapitulates the observations from Mostafavi et al. Nat Genet 2023 with regards to
#' the distance between variants and genes of cis-eQTLs.
#' The code here extends this to also calculate this for trans-eQTLs

# Load libraries
library(tidyverse)
library(data.table)
library(arrow)
library(extrafont)
library(tidyr)
library(broom)
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

# Declare function definitions

plot_to_pdf <- function(plot, path, width = 6, height = 6) {
  pdf(path,
      width = width, height = height, useDingbats = F)
  par(xpd = NA)

  print(plot)

  dev.off()
  embed_fonts(path, outfile=path)
}


# Define a function to compare case (Actual) vs control SNPs for each class using the 10kb cut-off
compare_case_control_10kb_0ld <- function(data) {
  results <- data %>%
    group_by(class) %>%
    summarize(
      actual_within_10kb = sum(within_10kb[control == "Actual"], na.rm = TRUE),   # Count of Actual SNPs within 10kb
      actual_total = sum(control == "Actual", na.rm = TRUE),                     # Total Actual SNPs
      control_within_10kb = sum(within_10kb[control == "Control"], na.rm = TRUE), # Count of Control SNPs within 10kb
      control_total = sum(control == "Control", na.rm = TRUE),                   # Total Control SNPs
      fisher_p_value = fisher.test(matrix(c(actual_within_10kb, actual_total - actual_within_10kb,
                                            control_within_10kb, control_total - control_within_10kb),
                                          nrow = 2))$p.value  # Fisher's Exact Test
    )
  return(results)
}

# Function to create a 2x2 contingency table
create_contingency_table <- function(actual_within_10kb, actual_total, control_within_10kb, control_total) {
  print(paste(actual_within_10kb, actual_total, control_within_10kb, control_total))
  contigency_matrix <- matrix(c(actual_within_10kb, actual_total - actual_within_10kb,
           control_within_10kb, control_total - control_within_10kb),
         nrow = 2)
  print(contigency_matrix)
  return(contigency_matrix)
}

# Function to perform Fisher's exact test and extract odds ratio and p-value
perform_fisher_test <- function(contingency_table) {
  fisher_result <- fisher.test(contingency_table)
  fisher_out <- list(odds_ratio = fisher_result$estimate, p_value = fisher_result$p.value)
  print(fisher_out)
  return(fisher_out)
}

# Function to calculate the standard error of the odds ratio
calculate_se_or <- function(actual_within_10kb, actual_total, control_within_10kb, control_total, odds_ratio) {
  # Standard error of the log odds ratio
  se_log_or <- sqrt(1 / actual_within_10kb + 1 / (actual_total - actual_within_10kb) +
                      1 / control_within_10kb + 1 / (control_total - control_within_10kb))

  # Standard error of the odds ratio
  se_or <- odds_ratio * se_log_or
  return(list("se_or" = se_or, "se_log_or" = se_log_or))
}

# Main function to compare case vs control SNPs within 10kb for each class
compare_case_control_10kb <- function(data) {

  results <- data %>%
    group_by(class) %>%
    summarize(
      actual_within_10kb = sum(within_10kb[control == "Actual"], na.rm = TRUE),   # Count of Actual SNPs within 10kb
      actual_total = sum(!within_10kb[control == "Actual"], na.rm = TRUE),                     # Total Actual SNPs
      control_within_10kb = sum(within_10kb[control == "Control"], na.rm = TRUE), # Count of Control SNPs within 10kb
      control_total = sum(!within_10kb[control == "Control"], na.rm = TRUE)                   # Total Control SNPs
    ) %>%
    rowwise() %>%
    mutate(
      # Create the contingency table
      contingency_table = list(create_contingency_table(actual_within_10kb, actual_total, control_within_10kb, control_total)),

      # Perform Fisher's exact test and extract odds ratio and p-value
      fisher_results = list(perform_fisher_test(contingency_table)),
      odds_ratio = fisher_results$odds_ratio,
      fisher_p_value = fisher_results$p_value,

      # Calculate standard error of the odds ratio
      se_or = list(calculate_se_or(actual_within_10kb, actual_total, control_within_10kb, control_total, odds_ratio))
    ) %>%
    unnest_wider(se_or) %>%
    ungroup()

  return(results)
}


# Example usage:
# Assuming `your_data` is the data frame you're working with
# results <- compare_case_control_10kb(your_data)
# print(results)

# Define a function to compare case (Actual) vs control SNPs for each class
compare_case_control <- function(data) {
  results <- data %>%
    group_by(class) %>% # Group by class (e.g., "GWAS hits", "cis-eQTL", etc.)
    mutate(abs_d_tss = abs(d_tss)) %>%
    summarize(
      actual_mean = mean(abs_d_tss[control == "Actual"], na.rm = TRUE),     # Mean distance for actual SNPs
      control_mean = mean(abs_d_tss[control == "Control"], na.rm = TRUE),   # Mean distance for control SNPs
      wilcox_p_value = wilcox.test(abs_d_tss[control == "Actual"],           # Wilcoxon test for case vs control
                                   abs_d_tss[control == "Control"])$p.value
    )
  return(results)
}

# Function to lift genome build to hg38 using rtracklayer (LiftOver)
lift_to_build_38 <- function(gwas_data, liftover_chain_file) {
  # Create a GRanges object for liftover
  gwas_data[, index := seq_along(gwas_data$SNP)]

  gr <- GRanges(
    seqnames = Rle(paste0("chr", str_remove(gwas_data$chromosome, "^chr"))),
    ranges = IRanges(start = gwas_data$bp, end = gwas_data$bp),
    strand = Rle("*"),
    index = gwas_data$index
  )

  # Perform liftover
  chain <- import.chain(liftover_chain_file)
  lifted_gr <- liftOver(gr, chain)

  # Flatten the result back to the data.table
  new_coords <- unlist(lifted_gr)
  new_coords_dt <- as.data.table(new_coords) %>% group_by(index) %>%
    filter(n() == 1) %>% as.data.table()

  gwas_data_filtered <- gwas_data[new_coords_dt[["index"]], , ]
  gwas_data_filtered[, chromosome := sub("chr", "", new_coords_dt[["seqnames"]])]
  gwas_data_filtered[, bp := new_coords_dt[["start"]]]
  gwas_data_filtered[ , index := NULL, ]

  return(gwas_data_filtered)
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

  phase_2_lead_eqtls_path <- variables["independent_variant_set"]
  results_phase2 <- fread(phase_2_lead_eqtls_path) %>%
    mutate(Z = beta/standard_error,
           SNP_MAPPING_ID = paste(chromosome, bp, pmin(eff_allele, non_eff_allele), pmax(eff_allele, non_eff_allele), sep="_")) %>%
    group_by(SNP_MAPPING_ID) %>%
    filter(n_distinct(variant) == 1)

  variant_reference <- read_parquet(variables["variants_with_unique_mapping_ids"])

  # external data
  base_external <- variables["mostafavi_zenodo"]

  snp_annotations <- fread(file.path(base_external, "snp_annotations/filter_snps.txt")) %>%
    separate_wider_delim(SNP, names=c("chromosome", "bp", "a1", "a2"), delim=":", cols_remove = F) %>%
    mutate(bp = as.integer(bp))

  snp_annotations_lifted <- lift_to_build_38(
    as.data.table(snp_annotations),
    liftover_chain_file = variables["liftover_to_b38"])

  snp_annotations_proc <- as_tibble(snp_annotations_lifted) %>%
    mutate(SNP_MAPPING_ID = paste(chromosome, bp, pmin(a1, a2), pmax(a1, a2), sep="_")) %>%
    filter(SNP_MAPPING_ID %in% variant_reference$mapping_id)

  gene_annotations <- fread(file.path(base_external, "gene_annots/pc_genes.txt"))
  gwas_hits <- fread(file.path(base_external, "gwas_hits.all_phenos.clumped_sorted.txt"), col.names = c("trait", "SNP", "p"))

  # Perform method


  # For every variant

  # Confounders:
  # Annotate eQTL/GWAS variants with MAF
  # Annotate eQTL/GWAS variants with LD score
  # Annotate eQTL/GWAS variants with gene density

  # Outcomes:
  # Annotate every eQTL/GWAS variant with distance to nearest TSS
  # Annotate every eQTL/GWAS gene with PLI

  # Optional others:
  # Regulation complexity
  # Functional annotations
  # Regulatory annotations

  # SNP selection:
  # biallelic autosomal SNPs, MAF > 0.01 in relevant ancestry
  # polymorphic in 1000G phase 3 (based on A. Prices LD score regerssion data)
  # Removed SNPs in LD to missense SNPs or protein truncating SNPs.
  # No MHC (chr6:28477797-33448354)

  eqtl_hits_annotated <- results_phase2 %>%
    inner_join(snp_annotations_proc, by=c("SNP_MAPPING_ID" = "SNP_MAPPING_ID"), relationship = "many-to-many")

  gwas_hits_annotated <- gwas_hits %>%
    inner_join(snp_annotations_proc, by=c("SNP" = "SNP"), relationship = "many-to-one")

  combined_results <- bind_rows(list("eQTL" = eqtl_hits_annotated, "GWAS" = gwas_hits_annotated), .id = "superclass") %>%
    mutate(class = case_when(
      superclass == "GWAS" ~ "GWAS hits",
      superclass == "eQTL" & type == "cis" ~ "cis-eQTL",
      superclass == "eQTL" & type == "trans" ~ "trans-eQTL")) %>%
    select(c("SNP", "class", "MAF", "L2", "TSSD", "d_tss"))
  # From every group, select eQTLs/GWAS hits with approximate equal MAF, LD score and gene density

  max_maf_delta <- 0.02
  max_ldsc_delta <- sd(snp_annotations$L2, na.rm=T)*0.1

  # Convert tibbles to data.tables for efficiency
  setDT(combined_results)
  setDT(snp_annotations_proc)

  combined_results[, `:=`(
    MAF_lower = MAF - max_maf_delta,
    MAF_upper = MAF + max_maf_delta,
    L2_lower = L2 - max_ldsc_delta,
    L2_upper = L2 + max_ldsc_delta
  )]

  setkey(snp_annotations_proc, MAF, L2, TSSD)
  setkey(combined_results, MAF_lower, MAF_upper, L2_lower, L2_upper, TSSD)

  # Join the datasets based on the precomputed ranges
  matched_snps <- combined_results[snp_annotations_proc, on = .(
    MAF_lower <= MAF, MAF_upper >= MAF,           # MAF within range
    L2_lower <= L2, L2_upper >= L2,               # LD score within range
    TSSD == TSSD                                  # Gene density (TSSD) exact match
  ), nomatch = 0L, allow.cartesian = TRUE]

  # Group by each focal SNP (SNP from combined_results) and sample 1,000 SNPs per group with replacement
  bootstrapped_results <- matched_snps[, .SD[sample(.N, 1000, replace = T)], by = .(SNP, class)]

  saveRDS(bootstrapped_results, "control_snps.rds")
  bootstrapped_results <- readRDS("control_snps.rds")

  control_snps <- bootstrapped_results %>%
    mutate(MATCHED_SNP = SNP_MAPPING_ID) %>%
    filter(i.SNP_MAPPING_ID != MATCHED_SNP) %>%
    select(SNP_MAPPING_ID = i.SNP_MAPPING_ID, SNP = i.SNP, class, MATCHED_SNP)

  #
  # control_snps <- fread("~/Downloads/control_snps.txt.gz")
  # control_snps <- fread("control_snps.txt.gz")

  control_snps_merged <- merge(snp_annotations_proc, control_snps, by = c("SNP_MAPPING_ID","SNP"))

  fwrite(control_snps_merged, "control_snps.txt.gz", sep="\t", col.names=T, row.names=F, quote=F)
  control_snps_merged <- fread("control_snps.txt.gz")

  combined_results_with_controls <- bind_rows(
    list("Actual" = combined_results,
         "Control" = control_snps_merged), .id="control") %>%
    mutate(
      class = if_else(is.na(class), "trans-eQTL", class),
      class_nested = case_when(control == "Actual" ~ class,
                             control == "Control" ~ paste0(class, " - matched SNPs")))

  combined_results_with_controls$abs_d_tss <- abs(combined_results_with_controls$d_tss)
  combined_results_with_controls$within_10kb <- combined_results_with_controls$abs_d_tss <= 10000
  comparison_10kb_new <- compare_case_control_10kb(combined_results_with_controls)

  delta_log_or <- (
    log(comparison_10kb_new[comparison_10kb_new$class == "cis-eQTL", "odds_ratio"])
    - log(comparison_10kb_new[comparison_10kb_new$class == "trans-eQTL", "odds_ratio"]))

  delta_log_or_se <- sqrt(
    comparison_10kb_new[comparison_10kb_new$class == "cis-eQTL", "se_log_or"]^2
      + comparison_10kb_new[comparison_10kb_new$class == "trans-eQTL", "se_log_or"]^2)

  delta_log_or_z <- delta_log_or / delta_log_or_se

  pnorm(as.numeric(delta_log_or_z), lower.tail=F)


  cis_trans_color_scale <- c("cis" = "#CC79A7", "trans" = "#009E73")

  color_manual <- c("GWAS hits" = "#56B4E9",
                    "GWAS hits - matched SNPs" = "#B1E0F9",
                    "trans-eQTL" = "#009E73",
                    "trans-eQTL - matched SNPs" = "#9FE8D4",
                    "cis-eQTL" = "#CC79A7",
                    "cis-eQTL - matched SNPs" = "#F1D1E3")

  # Most basic error bar
  ors <- ggplot(comparison_10kb_new) +
    geom_bar(aes(x=class, y=odds_ratio, fill = class), stat="identity") +
    scale_fill_manual(values = color_manual) +
    geom_errorbar(aes(x=class, ymin=odds_ratio - se_or, ymax=odds_ratio + se_or), width=0.4, colour="black", alpha=0.9, linewidth=1.3) +
    theme(legend.position="none", axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x = element_blank())

  ggsave(file.path(".", sprintf(
    "mostafavi_bar_chart.%s.pdf", format(Sys.Date(), "%Y%m%d"))), ors, width=80, height=60, units = "mm")
  plot_to_pdf(
    ors, file.path(".", sprintf(
      "mostafavi_bar_chart.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 3, height = 7)

  # Plot

  # distance to nearest TSS

  # Distance to nearest gene
  bin_definitions <- seq(-1e5, +1e5, by= 1e4)
  midpoints <- (bin_definitions[-length(bin_definitions)] + bin_definitions[-1]) / 2

  tss_distances <- combined_results_with_controls %>%
    filter(!is.na(class)) %>%
    ungroup() %>%
    mutate(distance_bin = as.character(cut(d_tss, breaks = bin_definitions, include.lowest = T, right = F, labels=midpoints))) %>%
    group_by(class, class_nested, control) %>%
    mutate(total_snps = n()) %>%
    ungroup() %>%
    group_by(distance_bin, class, class_nested, control) %>%
    summarise(prop_of_snps = n()/dplyr::first(total_snps)) %>%
    ungroup() %>%
    mutate(distance_bin = as.numeric(distance_bin))

  distance_plot <- ggplot(
    tss_distances,
    aes(x=distance_bin, y=prop_of_snps, colour=class_nested)) +
    geom_vline(xintercept = 0, linewidth = (0.5 / (ggplot2::.pt * 72.27/96)), colour="grey80") +
    geom_line(linewidth=(0.5 / (ggplot2::.pt * 72.27/96))) +
    geom_point(size=(1 / (ggplot2::.pt * 72.27/96)), shape=16) +
    scale_colour_manual(values = color_manual) +
    scale_x_continuous("Distance to nearest TSS", breaks=c(-1e5, -0.5e5, 0, 0.5e5, 1e5), labels=paste0(c(-100, -50, 0, 50, 100), "kb"))+
    scale_y_continuous("Proportion of SNPs") +
    facet_grid(cols=vars(class)) +
    guides(
      colour = guide_legend(position = "inside"),
    ) +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

  distance_plot


  ggsave(file.path(".", sprintf(
    "mostafavi_tss.%s.pdf", format(Sys.Date(), "%Y%m%d"))), distance_plot, width=80, height=50, units = "mm")

  plot_to_pdf(
    distance_plot, file.path(".", sprintf(
      "mostafavi_tss.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 10, height = 8)


  combined_results %>%
    group_by(type) %>%
    summarise(median_pli = median(pLI, na.rm=T),
              mean_pli = mean(pLI, na.rm=T),
              proportion_of_high_pli_genes = sum(pLI > 0.9, na.rm = T)/n())

  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}