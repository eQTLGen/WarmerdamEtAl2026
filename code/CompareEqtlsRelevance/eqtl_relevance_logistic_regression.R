#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)
library(arrow)


# Declare constants
old <- theme_set(theme_classic(base_size = 8))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line("grey90"),
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
  featureDescriptions <- fread("eQTLRelevanceFeatureDescriptions.txt")

  # Variant reference
  variant_reference <- read_parquet(variables["variant_reference"])
  variant_index_to_variant_id <- variant_reference %>% select(variant_index, variant) %>% deframe()

  # Indepedent variant set
  independent_variant_set_annot_filtered <- fread(variables["independent_variant_set"])

  # Number of cis-eQTL credible sets
  additional_annots <- independent_variant_set_annot_filtered %>%
    filter(!is.na(SusieRss_CS), type == "cis") %>%
    group_by(phenotype) %>%
    mutate(GeneNCisCredibleSets = n(), beta_abs = abs(beta)) %>%
    select(variant, phenotype, n_unique_CSs, beta_abs)


  # Process input
  ## Get features

  renaming_map <- c("SNPATACseqOpenness" = "SNPATACseqMultiome")
  renaming_map2 <- names(renaming_map)
  names(renaming_map2) <- unname(renaming_map)

  subset_significant_features <- fread("testdata2.txt.SubsetSignificantFeatures.txt") %>%
    separate_wider_delim(cols = c(`-`), delim = "-", names = c("variant", "phenotype")) %>%
    rename_with(.fn = ~ str_replace_all(.x, "-", "_")) %>%
    rename(renaming_map2) %>%
    inner_join(additional_annots)

  fwrite(subset_significant_features %>% select(-eQTLscEQTLGenBetaStdev, -beta_abs), "testdata2.txt.SubsetSignificantFeatures5.txt", sep="\t", quote=F, col.names=T, row.names=F)
  fwrite(subset_significant_features %>% select(-eQTLscEQTLGenBetaStdev, -eQTLscMetaBrainBetaStdev, -beta_abs), "testdata_20260724_SubsetSignificantFeatures_noScMetaBrain.txt", sep="\t", quote=F, col.names=T, row.names=F)

  gene_matched_feature_table <- fread("FunctionallyActiveVsIrrelevantCisEQTLs.txt.SubsetSignificantFeatures.txt") %>%
    separate_wider_delim(cols = c(`-`), delim = "-", names = c("variant", "phenotype")) %>%
    rename_with(.fn = ~ str_replace_all(.x, "-", "_")) %>%
    rename(renaming_map2) %>%
    inner_join(additional_annots)


  any_na <- subset_significant_features %>% summarise(across(everything(), ~ any(is.infinite(.x)))) %>%
    as.vector()

  print(any_na)

  any_na_gene_matched <- gene_matched_feature_table %>% summarise(across(everything(), ~ any(is.infinite(.x)))) %>%
    as.vector()

  print(any_na_gene_matched)

  all_features <- names(subset_significant_features)

  predictors <- setdiff(
    all_features,
    c("variant", "phenotype", "Relevance", "eQTLscEQTLGenBetaStdev", "beta_abs")
  )

  correlation_matrix <- cor(subset_significant_features[,c(predictors, "Relevance")], method="spearman")
  correlation_df <- as.data.table(correlation_matrix, keep.rownames = "Feature") %>%
    pivot_longer(-Feature, names_to="Feature2", values_to = "Spearman")


  # Cluster features using hierarchical clustering
  feature_order <- hclust(
    as.dist(1 - correlation_matrix),
    method = "complete"
  )$order

  clustered_features <- rownames(correlation_matrix)[feature_order]

  # Apply clustered ordering
  correlation_df <- correlation_df %>%
    mutate(
      Feature  = factor(Feature, levels = clustered_features),
      Feature2 = factor(Feature2, levels = rev(clustered_features))
    )

  # Plot heatmap
  ggplot(correlation_df, aes(x = Feature2, y = Feature, fill = Spearman)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Spearman\nrho"
    ) +
    coord_fixed() +
    labs(
      x = NULL,
      y = NULL) +
    theme_minimal(base_size = 6) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )

  ggsave("correlation_plot.pdf", width=180, height=200, units="mm")

  # Assess effect of CS count on Borzoi association

  predictors_no_cs_count <- setdiff(
    all_features,
    c("variant", "phenotype", "Relevance", "eQTLscEQTLGenBetaStdev", "beta_abs", "n_unique_CSs")
  )

  relevance_formula_no_cs_count <- as.formula(
    paste("Relevance ~", paste(predictors_no_cs_count, collapse = " + "))
  )

  relevance_mod_no_cs_count <- glm(
    relevance_formula_no_cs_count,
    data = subset_significant_features,
    family = binomial(link = "logit")
  )

  summary(relevance_mod_no_cs_count)


  relevance_formula <- as.formula(
    paste("Relevance ~", paste(predictors, collapse = " + "))
  )

  full_relevance_mod <- glm(
    relevance_formula,
    data = subset_significant_features,
    family = binomial(link = "logit")
  )

  summary(full_relevance_mod)
  res_full <- as_tibble(summary(full_relevance_mod)$coefficients, rownames="Feature")

  gene_matched_relevance_mod <- glm(
    relevance_formula,
    data = gene_matched_feature_table,
    family = binomial(link = "logit")
  )

  summary(gene_matched_relevance_mod)
  res_gene_matched <- as_tibble(summary(gene_matched_relevance_mod)$coefficients, rownames="Feature")

  full_relevance_mod_stdized <- glm(
    relevance_formula,
    data = subset_significant_features %>% mutate(across(all_of(predictors), .fn = scale)),
    family = binomial(link = "logit")
  )

  summary(full_relevance_mod_stdized)
  res_full_stdized <- as_tibble(summary(full_relevance_mod_stdized)$coefficients, rownames="Feature")

  full_results_relevance <- bind_rows(
    full = res_full,
    full_std = res_full_stdized,
    gene_matched = res_gene_matched, .id="Gene") %>%
    left_join(featureDescriptions, by =c("Feature"="Feature"))

  fwrite(full_results_relevance, "logistic_regression_model_eqtl_relevance.txt.gz",
         sep="\t", quote=F, row.names=F, col.names=T)

  res_full_stdized_proc <- res_full_stdized %>%
    filter(Feature != "(Intercept)" & `Pr(>|z|)` < 0.05) %>%
    arrange(abs(Estimate)) %>%
    left_join(featureDescriptions) %>%
    mutate(
      Feature2 = sprintf("%s:%s", Level, str_remove(Feature, "eQTL|SNP|Gene")),
      Feature2 = ordered(Feature2, levels = Feature2))

  max_abs <- max(abs(res_full_stdized_proc$Estimate) + res_full_stdized_proc$`Std. Error`)

  ggplot(res_full_stdized_proc) +
    geom_bar(aes(x=Feature2, y=Estimate, fill=if_else(sign(Estimate)==1, "P", "N")), stat="identity") +
    geom_errorbar( aes(x=Feature2, ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`), width=0.4, colour="black", size=2/3/.pt) +
    coord_flip() +
    scale_y_continuous(name="Standardized beta coefficient", breaks=seq(-0.5, 0.5, 0.1)) +
    guides(
      x = guide_axis(cap = "both"), # Cap both ends
      y = guide_axis(cap = "both") # Cap the upper end
    ) + theme(axis.title.y = element_blank(), legend.position = "none", axis.text=element_text(size=7))

  ggsave("actionable_eqtls_standardized_effect_sizes.pdf", width=92, height=104, units="mm")

  average_distance_to_gene <- subset_significant_features %>% group_by(Relevance) %>%
    filter(SNPGeneDistance != 0) %>%
    mutate(relevance_label = if_else(Relevance==1,
                                     "Cis-eQTL with trans-eQTLs",
                                     "Cis-eQTL without trans-eQTLs")) %>%
    summarise(medianSnpGeneDistance = median(SNPGeneDistance), meanSnpGeneDistance=mean(SNPGeneDistance))

  average_distance_to_gene_2 <- subset_significant_features %>% group_by(Relevance) %>%
    mutate(relevance_label = if_else(Relevance==1,
                                     "Cis-eQTL with trans-eQTLs",
                                     "Cis-eQTL without trans-eQTLs")) %>%
    summarise(medianSnpGeneDistance = median(SNPGeneDistance), meanSnpGeneDistance=mean(SNPGeneDistance))

  gene_matched_feature_table %>% group_by(Relevance) %>%
    filter(SNPGeneDistance != 0) %>%
    mutate(relevance_label = if_else(Relevance==1,
                                     "Cis-eQTL with trans-eQTLs",
                                     "Cis-eQTL without trans-eQTLs")) %>%
    summarise(medianSnpGeneDistance = median(SNPGeneDistance), meanSnpGeneDistance=mean(SNPGeneDistance))


  t.test(gene_matched_feature_table$SNPGeneDistance[gene_matched_feature_table$SNPGeneDistance != 0 & gene_matched_feature_table$Relevance==1],
         gene_matched_feature_table$SNPGeneDistance[gene_matched_feature_table$SNPGeneDistance != 0 & gene_matched_feature_table$Relevance==0])


  wilcox.test(
    subset_significant_features$SNPGeneDistance[subset_significant_features$SNPGeneDistance != 0 & subset_significant_features$Relevance==1],
    subset_significant_features$SNPGeneDistance[subset_significant_features$SNPGeneDistance != 0 & subset_significant_features$Relevance==0])

  wilcox.test(
    subset_significant_features$SNPGeneDistance[subset_significant_features$Relevance==1],
    subset_significant_features$SNPGeneDistance[subset_significant_features$Relevance==0])

  t.test(subset_significant_features$SNPGeneDistance[subset_significant_features$SNPGeneDistance != 0 & subset_significant_features$Relevance==1],
         subset_significant_features$SNPGeneDistance[subset_significant_features$SNPGeneDistance != 0 & subset_significant_features$Relevance==0])

  mbreaks_non_null <- seq(0,1e6, 5e4)
  mbreaks <- c(-Inf, mbreaks_non_null)
  breaks_formatted <- sprintf("%.0fkb", mbreaks_non_null / 1000)
  breaks_formatted[length(breaks_formatted)] <- "1Mb"
  lables <- c("0", paste(lag(breaks_formatted), breaks_formatted, sep="-")[2:length(breaks_formatted)])
  bin_midpoints <- (mbreaks_non_null + lag(mbreaks_non_null)) / 2
  bin_midpoints[1] <- 0
  names(bin_midpoints) <- lables

full_feature_table_breaks <- subset_significant_features %>%
    select(Relevance, SNPGeneDistance) %>%
    mutate(relevance_label = if_else(Relevance==1,
                                     "Cis-eQTL with trans-eQTLs",
                                     "Cis-eQTL without trans-eQTLs"),
           bins = cut(SNPGeneDistance, breaks = mbreaks, ordered_result=T, labels=lables),
           bin_midpoint = bin_midpoints[bins]
    ) %>%
    group_by(bins, bin_midpoint, Relevance, relevance_label) %>%
    summarise(n_in_group_bin_combo = n()) %>%
    group_by(Relevance) %>%
    mutate(prop = n_in_group_bin_combo / sum(n_in_group_bin_combo)) %>%
    ungroup()


  gene_matched_feature_table_breaks <- gene_matched_feature_table %>%
    select(Relevance, SNPGeneDistance) %>%
    mutate(relevance_label = if_else(Relevance==1,
                                     "Cis-eQTL with trans-eQTLs",
                                     "Cis-eQTL without trans-eQTLs"),
           bins = cut(SNPGeneDistance, breaks = mbreaks, ordered_result=T, labels=lables),
           bin_midpoint = bin_midpoints[bins]
    ) %>%
    group_by(bins, bin_midpoint, Relevance, relevance_label) %>%
    summarise(n_in_group_bin_combo = n()) %>%
    group_by(Relevance) %>%
    mutate(prop = n_in_group_bin_combo / sum(n_in_group_bin_combo)) %>%
    ungroup()


  ggplot(full_feature_table_breaks, aes(x=bin_midpoint,y=prop, color=relevance_label)) +
    geom_vline(data = average_distance_to_gene, aes(xintercept = meanSnpGeneDistance), linewidth=(0.5 / (ggplot2::.pt * 72.27/96))) +
    geom_point() + geom_line(aes(group = relevance_label)) +
    theme(legend.position = "null", panel.grid.major=element_line(color="grey90"), axis.text=element_text(size=7)) +
    scale_y_continuous(name = "Proportion of variants", limits=c(0,0.4), breaks=seq(0,0.4,0.05), labels=sprintf("%.0f%%", seq(0, 40, 5))) +
    scale_x_continuous(limits=c(-1,1e6), breaks=c(0,2e5,4e5,6e5,8e5,1e6), labels =sprintf("%.0f kbp", c(0,2e5,4e5,6e5,8e5,1e6)/1000)) +
    guides(
      x = guide_axis(cap = "both"), # Cap both ends
      y = guide_axis(cap = "both") # Cap the upper end
    )

  ggsave("eqtl_gene_distance_effect_sizes.pdf", width=74, height=66, units="mm")

  ggplot(gene_matched_feature_table_breaks, aes(x=bin_midpoint,y=prop, color=relevance_label)) +
    geom_point() + geom_line(aes(group = relevance_label)) + theme(legend.position = "bottom", axis.text.x=element_text(angle=45, hjust=1)) +
    scale_y_log10(name = "Proportion of variants")

  # Do promotors vs enhancers more often overlap cis-eQTLs with trans-eQTLs

  fisher.test(
    table(subset_significant_features[c("Relevance", "SNPENCODE_dELS")]))

  minimal_predictors <- c(
    "SNPMAF",
    "eQTLAbsZ",
    "eQTLFinemapPIP",
    "eQTLCredibleSetRank",
    "SNPLDScore",
    "eQTLFinemappingLambda",
    "eQTLHeterogeneityI2",
    "GeneDensity",
    "GeneGCContent",
    "SNP5mer",
    "GeneNrTranscripts",
    "SNPisIndel",
    "n_unique_CSs"
  )

  relevance_formula_enhancers_promoters <- as.formula(
    paste("Relevance ~", paste(c(minimal_predictors, "SNPENCODE_PLS", "SNPENCODE_pELS", "SNPENCODE_dELS"), collapse = " + "))
  )

  enhancers_promoters_mod <- glm(
    relevance_formula_enhancers_promoters,
    data = subset_significant_features,
    family = binomial(link = "logit")
  )

  summary(enhancers_promoters_mod)
  enhancers_promotors_subset <- as_tibble(summary(enhancers_promoters_mod)$coefficients, rownames="Feature")
  fwrite(enhancers_promotors_subset, "enhancers_promotors_logistic_regression_results.txt.gz", sep="\t", col.names=T, row.names=F, quote=F)

  # Overlap relevant / irrelevant cis-eQTLs with experimentally validated eQTLs

  experimental_data_4 <- fread('SupplementaryData 2/Table_S4.txt')

  variants <- c(unlist(str_split(experimental_data_4$SCANDAL_locus, "_")), subset_significant_features$variant)
  fwrite(tibble(vars=variants), "all_variants_to_query_for_ld.txt", col.names=F, row.names=F)

  ld_matrices <- readRDS("tmp_extract_ld/ld_matrices.rds")

  lookup <- map_dfr(names(ld_matrices), function(chr) {
    cat(chr)

    ld <- ld_matrices[[chr]]^2
    rownames(ld) <- variant_index_to_variant_id[rownames(ld)]
    colnames(ld) <- variant_index_to_variant_id[colnames(ld)]

    query <- intersect(
      subset_significant_features$variant,
      rownames(ld)
    )

    tested <- intersect(
      experimental_data_4$SCANDAL_locus,
      colnames(ld)
    )

    ld_sub <- ld[query, tested, drop = FALSE]

    str(ld_sub)

    idx <- apply(ld_sub, 1, which.max)

    tibble(
      variant = query,
      paired_variant = tested[idx],
      R2 = ld_sub[cbind(seq_along(query), idx)]
    )
  })

  lookup_2 <- lookup %>%
    mutate(
      variant_to_pair = case_when(
        variant %in% experimental_data_4$SCANDAL_locus ~ variant,
        R2 >= 0.8 ~ paired_variant,
        TRUE ~ NA_character_
      ),
      R2 = if_else(
        variant %in% experimental_data_4$SCANDAL_locus,
        2,
        R2
      )
    )

  subset_significant_features_with_b_cell_crispr <- subset_significant_features %>%
    left_join(lookup_2, by = "variant") %>%
    left_join(
      experimental_data_4,
      by = c(
        "variant_to_pair" = "SCANDAL_locus",
        "phenotype" = "Target_Ensembl"
      )
    ) %>%
    mutate(experimantally_validated = if_else(!is.na(p_adj_BH), 1, 0))

  fwrite(subset_significant_features_with_b_cell_crispr, "features-with_b_cell_crispr.txt.gz",
         quote=F, sep="\t", row.names=F, col.names=T)

  confusion_table_1 <- subset_significant_features_with_b_cell_crispr %>%
    group_by(Relevance, experimantally_validated) %>% tally()
  confusion_table_1
  confusion_table_2 <- confusion_table_1 %>%
    pivot_wider(values_from = n, names_from = experimantally_validated) %>%
    as.data.table() %>%
    as.matrix(rownames="Relevance")

  fisher.test(confusion_table_2)

  confusion_table_3 <- confusion_table_2
  colnames(confusion_table_3) <- c("No", "Yes")
  rownames(confusion_table_3) <- c("No", "Yes")


  library(caret)
  confusionMatrix(confusion_table_3, positive = "Yes")

  predictors_experimentally_validated <- setdiff(
    all_features,
    c("variant", "phenotype", "eQTLscEQTLGenBetaStdev", "beta_abs")
  )

  predictors_experimentally_validated_2 <- c(minimal_predictors, "Relevance")


  validated_formula <- as.formula(
    paste("experimantally_validated ~", paste(predictors_experimentally_validated, collapse = " + "))
  )

  validated_formula_2 <- as.formula(
    paste("experimantally_validated ~", paste(predictors_experimentally_validated_2, collapse = " + "))
  )

  validated_mod <- glm(
    validated_formula,
    data = subset_significant_features_with_b_cell_crispr,
    family = binomial(link = "logit")
  )

  summary(validated_mod)


  validated_mod_2 <- glm(
    validated_formula_2,
    data = subset_significant_features_with_b_cell_crispr,
    family = binomial(link = "logit")
  )

  summary(validated_mod_2)



  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}