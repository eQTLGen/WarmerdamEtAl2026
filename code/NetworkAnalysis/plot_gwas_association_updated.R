#!/usr/bin/env Rscript
#' Generate figure 5a and b

# Load libraries
library(tidyverse)
library(data.table)
library(ggrepel)
library(patchwork)

# Declare constants
old <- theme_set(theme_classic(base_size = 8))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  panel.grid.minor = element_blank(),
  text = element_text(family="Helvetica", colour="black"),
  title = element_text(colour = "#595A5C", face = "plain")
)

trait_categories <- c(
  # Immune mediated diseases
  "PASS_CD_deLange2017" = "Immune mediated diseases",
  "PASS_Celiac" = "Immune mediated diseases",
  "PASS_IBD_deLange2017" = "Immune mediated diseases",
  "PASS_Lupus" = "Immune mediated diseases",
  "PASS_Lupus_Wang2021" = "Immune mediated diseases",
  "PASS_Multiple_sclerosis2019" = "Immune mediated diseases",
  "PASS_Primary_biliary_cirrhosis" = "Immune mediated diseases",
  "PASS_Rheumatoid_Arthritis" = "Immune mediated diseases",
  "PASS_UC_deLange2017" = "Immune mediated diseases",
  "UKB_460K.disease_AID_ALL" = "Immune mediated diseases",
  "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED" = "Immune mediated diseases",
  "UKB_460K.disease_ASTHMA_DIAGNOSED" = "Immune mediated diseases",
  "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP" = "Immune mediated diseases",
  "PASS_Type_1_Diabetes" = "Immune mediated diseases",

  # Cell type composition
  "UKB_460K.blood_EOSINOPHIL_COUNT" = "Cell type composition",
  "UKB_460K.blood_LYMPHOCYTE_COUNT" = "Cell type composition",
  "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN" = "Cell type composition",
  "UKB_460K.blood_MONOCYTE_COUNT" = "Cell type composition",
  "UKB_460K.blood_PLATELET_COUNT" = "Cell type composition",
  "UKB_460K.blood_RBC_DISTRIB_WIDTH" = "Cell type composition",
  "UKB_460K.blood_RED_COUNT" = "Cell type composition",
  "UKB_460K.blood_WHITE_COUNT" = "Cell type composition",
  "Basophil_count" = "Cell type composition",
  "Basphil_percentage_of_leucocytes" = "Cell type composition",
  "Eosinophil_count" = "Cell type composition",
  "Eosinophil_percentage_of_white_cells" = "Cell type composition",
  "Hematocrit" = "Cell type composition",
  "Hemoglobin" = "Cell type composition",
  "White_blood_cell_count" = "Cell type composition",
  "Lymphocyte_percentage_of_white_cells" = "Cell type composition",
  "Mean_corpuscular_hemoglobin" = "Cell type composition",
  "Mean_corpuscular_hemoglobin_concentration" = "Cell type composition",
  "Mean_corpuscular_volume" = "Cell type composition",
  "Mean_reticulocyte_volume" = "Cell type composition",
  "Mean_spheric_corpuscular_volume" = "Cell type composition",
  "Monocyte_count" = "Cell type composition",
  "Neutrophil_count" = "Cell type composition",
  "Neutrophil_percentage_of_white_cells" = "Cell type composition",
  "Plateletcrit" = "Cell type composition",
  "Platelet_distribution_width" = "Cell type composition",
  "High_light_scatter_reticulocyte_count" = "Cell type composition",
  "High_light_scatter_reticulocyte_percentage_of_red_cells" = "Cell type composition",
  "Immature_fraction_of_reticulocytes" = "Cell type composition",
  "Monocyte_percentage_of_white_cells" = "Cell type composition",
  "Mean_platelet_volume" = "Cell type composition",

  # Metabolism & cardiometabolic disease
  "UKB_460K.impedance_BASAL_METABOLIC_RATEz" = "Metabolism & cardiometabolic disease",
  "PASS_FastingGlucose_Manning" = "Metabolism & cardiometabolic disease",
  "UKB_460K.biochemistry_Glucose" = "Metabolism & cardiometabolic disease",
  "UKB_460K.biochemistry_HbA1c" = "Metabolism & cardiometabolic disease",
  "PASS_Type_2_Diabetes" = "Metabolism & cardiometabolic disease",
  "PASS_AtrialFibrillation_Yuan2025" = "Metabolism & cardiometabolic disease",
  "PASS_Coronary_Artery_Disease" = "Metabolism & cardiometabolic disease",
  "UKB_460K.bp_DIASTOLICadjMEDz" = "Metabolism & cardiometabolic disease",
  "UKB_460K.bp_SYSTOLICadjMEDz" = "Metabolism & cardiometabolic disease",
  "UKB_460K.disease_CARDIOVASCULAR" = "Metabolism & cardiometabolic disease",
  "UKB_460K.disease_HYPERTENSION_DIAGNOSED" = "Metabolism & cardiometabolic disease",
  "UKB_460K.biochemistry_Cholesterol" = "Metabolism & cardiometabolic disease",
  "UKB_460K.biochemistry_HDLcholesterol" = "Metabolism & cardiometabolic disease",
  "UKB_460K.biochemistry_LDLdirect" = "Metabolism & cardiometabolic disease",
  "UKB_460K.biochemistry_Triglycerides" = "Metabolism & cardiometabolic disease",
  "High_Density_Lipoprotein_Cholesterol_levels" = "Metabolism & cardiometabolic disease",
  "Low_Density_Lipoprotein_Cholesterol_levels" = "Metabolism & cardiometabolic disease",
  "Non-HDL_Cholesterol_levels" = "Metabolism & cardiometabolic disease",
  "Total_Cholesterol_Measurement" = "Metabolism & cardiometabolic disease",
  "Triglyceride_Levels" = "Metabolism & cardiometabolic disease",

  # Neuropsychiatric & cognitive traits
  "PASS_ADHD_Demontis2023" = "Neuropsychiatric & cognitive traits",
  "PASS_Alzheimers_Bellenguez2022" = "Neuropsychiatric & cognitive traits",
  "PASS_BIP_Mullins2021" = "Neuropsychiatric & cognitive traits",
  "PASS_MDD_Als2023" = "Neuropsychiatric & cognitive traits",
  "PASS_Schizophrenia_Trubetskoy2022" = "Neuropsychiatric & cognitive traits",
  "UKB_460K.mental_NEUROTICISM" = "Neuropsychiatric & cognitive traits",
  "PASS_Insomnia_Watanabe2022" = "Neuropsychiatric & cognitive traits",
  "PASS_SleepDuration_Dashti2019" = "Neuropsychiatric & cognitive traits",
  "UKB_460K.other_MORNINGPERSON" = "Neuropsychiatric & cognitive traits",
  "PASS_Intelligence_SavageJansen2018" = "Neuropsychiatric & cognitive traits",
  "PASS_ReactionTime_Davies2018" = "Neuropsychiatric & cognitive traits",
  "PASS_VerbalNumericReasoning_Davies2018" = "Neuropsychiatric & cognitive traits",
  "PASS_Worry_Nagel2018" = "Neuropsychiatric & cognitive traits",
  "PASS_SWB" = "Neuropsychiatric & cognitive traits",
  "PASS_GeneralRiskTolerance_KarlssonLinner2019" = "Neuropsychiatric & cognitive traits",
  "PASS_DrinksPerWeek_Liu2019" = "Neuropsychiatric & cognitive traits",
  "UKB_460K.cov_SMOKING_STATUS" = "Neuropsychiatric & cognitive traits",

  # Other traits (anthropometric, reproductive, lifestyle, lung, cancer, pigmentation, education)
  "UKB_460K.body_BMIz" = "Other traits",
  "UKB_460K.body_HEIGHTz" = "Other traits",
  "UKB_460K.body_WHRadjBMIz" = "Other traits",
  "UKB_460K.body_BALDING1" = "Other traits",
  "UKB_460K.bmd_HEEL_TSCOREz" = "Other traits",
  "UKB_460K.repro_NumberChildrenEverBorn_Pooled" = "Other traits",
  "UKB_460K.repro_MENARCHE_AGE" = "Other traits",
  "UKB_460K.repro_MENOPAUSE_AGE" = "Other traits",
  "UKB_460K.pigment_HAIR" = "Other traits",
  "UKB_460K.cov_EDU_COLLEGE" = "Other traits",
  "UKB_460K.cov_EDU_YEARS" = "Other traits",
  "UKB_460K.cancer_BREAST" = "Other traits",
  "UKB_460K.disease_RESPIRATORY_ENT" = "Other traits",
  "UKB_460K.lung_FEV1FVCzSMOKE" = "Other traits",
  "UKB_460K.lung_FVCzSMOKE" = "Other traits"
)


disease_traits <- c(
  # Immune mediated diseases
  "PASS_CD_deLange2017",
  "PASS_Celiac",
  "PASS_IBD_deLange2017",
  "PASS_Lupus",
  "PASS_Lupus_Wang2021",
  "PASS_Multiple_sclerosis2019",
  "PASS_Primary_biliary_cirrhosis",
  "PASS_Rheumatoid_Arthritis",
  "PASS_UC_deLange2017",
  "UKB_460K.disease_AID_ALL",
  "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED",
  "UKB_460K.disease_ASTHMA_DIAGNOSED",
  "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP",
  "PASS_Type_1_Diabetes",

  # Metabolism & cardiometabolic disease
  "PASS_Type_2_Diabetes",
  "PASS_AtrialFibrillation_Yuan2025",
  "PASS_Coronary_Artery_Disease",
  "UKB_460K.disease_CARDIOVASCULAR",
  "UKB_460K.disease_HYPERTENSION_DIAGNOSED",

  # Neuropsychiatric & cognitive traits
  "PASS_ADHD_Demontis2023",
  "PASS_Alzheimers_Bellenguez2022",
  "PASS_BIP_Mullins2021",
  "PASS_MDD_Als2023",
  "PASS_Schizophrenia_Trubetskoy2022",

  "UKB_460K.cancer_BREAST",
  "UKB_460K.disease_RESPIRATORY_ENT"
)


# Declare function definitions

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

  gwas_finemap <- variables["gwas_signals"]
  coloc_glob <- variables["coloc_cleaned"]

  gwas_selection <- readxl::read_xlsx(variables["eqtlgen_gwass"]) %>%
    select(trait_identifier, trait_name)

  colocalizations <- bind_rows(mapply(fread, Sys.glob(coloc_glob), SIMPLIFY = F))

  colocalizations_with_lead_p <- colocalizations %>%
    mutate(eQTL_lead_eQTL_P = 2*pnorm(abs(eQTL_lead_eQTL_Z), lower.tail = F),
           eQTL_lead_GWAS_P = 2*pnorm(abs(eQTL_lead_GWAS_Z), lower.tail = F)
           )

  fwrite(colocalizations_with_lead_p, "eQTL_GWAS_coloc_cleaned_20260724.txt.gz", row.names=F, col.names=T, sep="\t", quote=F)

  colocalizations_proc <- colocalizations %>%
    filter(PP_H4 > 0.8) %>%
    distinct(GWAS, GWAS_locus, eQTL_type, signal2) %>%
    mutate(coloc = T) %>%
    pivot_wider(names_from = "eQTL_type", values_from=coloc, values_fill=F)

  joined <- fread(gwas_finemap) %>%
    mutate(gwas_signal = as.integer(str_extract(LBF_Column, "\\d+$"))) %>%
    separate_wider_delim(File, "__", names=c("GWAS", "GWAS_locus"), too_many = "drop") %>%
    mutate(GWAS_locus = str_replace(GWAS_locus, "_", "-")) %>%
    left_join(colocalizations_proc, by = c("GWAS" = "GWAS", "GWAS_locus" = "GWAS_locus", "gwas_signal" = "signal2")) %>%
    filter(!GWAS %in% c("PASS_Celiac", "PASS_Celiac_Alam2025"), Has_Signal_Lbf_P_lenient) %>%
    mutate(cis = replace_na(cis, F),
           trans = replace_na(trans, F))


  colocalization_stats_loci <- joined %>%
    group_by(GWAS, GWAS_locus) %>%
    summarise(cis = any(cis), trans = any(trans)) %>%
    group_by(GWAS) %>%
    summarise(
      n_signals = n(),
      n_coloc_cis = sum(cis),
      n_coloc_trans = sum(trans),
      n_coloc_all = sum(cis | trans)
    ) %>%
    pivot_longer(starts_with("n_coloc"),
                 names_pattern="n_coloc_(.+)",
                 values_to="n_coloc", names_to="eqtl_type") %>%
    mutate(
      trait_type = trait_categories[GWAS],
      disease_trait = GWAS %in% disease_traits,
      percentage = n_coloc / n_signals * 100,
      eqtl_type = ordered(eqtl_type, labels=c("cis-eQTLs", "trans-eQTLs", "all"), levels=c("cis", "trans", "all")))

  colocalization_stats <- joined %>%
    group_by(GWAS) %>%
    summarise(
      n_signals = n(),
      n_coloc_cis = sum(cis),
      n_coloc_trans = sum(trans),
      n_coloc_all = sum(cis | trans)
    ) %>%
    pivot_longer(starts_with("n_coloc"),
                 names_pattern="n_coloc_(.+)",
                 values_to="n_coloc", names_to="eqtl_type") %>%
    mutate(
      trait_type = trait_categories[GWAS],
      disease_trait = GWAS %in% disease_traits,
      percentage = n_coloc / n_signals * 100,
      eqtl_type = ordered(eqtl_type, labels=c("cis-eQTLs", "trans-eQTLs", "all"), levels=c("cis", "trans", "all"))) %>%
    inner_join(colocalization_stats_loci, by = c("GWAS", "trait_type", "eqtl_type", "disease_trait"), suffix=c("_cs", "_loci")) %>%
    left_join(gwas_selection, by = c("GWAS" = "trait_identifier")) %>%
    select(trait_name, everything())

    fwrite(colocalization_stats, "colocalization_stats.txt", sep="\t", col.names=T, row.names=F)

  trait_type_colors <- c(
    "Cell type composition" = "#D55E00",
    "Immune mediated diseases" = "#009E73",
    "Metabolism & cardiometabolic disease" = "#56B4E9",
    "Neuropsychiatric & cognitive traits" = "#CC79A7",
    "Other traits" = "#AAAAAA")

  p <- ggplot(colocalization_stats %>% filter(eqtl_type=="all"),
              aes(n_signals, n_coloc, color = trait_type, group=trait_type)) +
    geom_ribbon(stat='smooth', method = "lm", se=TRUE, alpha=0.1, linewidth=0, color="white") +
    geom_point(alpha=0.5, shape=16) +
    geom_text_repel(data = . %>% group_by(trait_type) %>% slice_max(n_signals), aes(label=trait_type), size=2) +
    geom_line(stat='smooth', method = "lm", alpha=1, linewidth=(0.5 / (ggplot2::.pt * 72.27/96))) +
    scale_colour_manual(values=trait_type_colors, name="Trait category", labels=set_names(str_wrap(names(trait_type_colors), 18), names(trait_type_colors))) +
    scale_x_continuous("Number of GWAS signals") +
    scale_y_continuous("Number of colocalizing GWAS signals") +
    guides(colour = guide_legend(ncol = 2))+
    theme(aspect.ratio = 1, legend.position = "none", legend.direction = "vertical")

  ggsave(sprintf("colocalizations_regression_%s.pdf", Sys.Date()), p, width=64, height=64, units="mm")
  #
  p2 <- ggplot(colocalization_stats, aes(eqtl_type, percentage, fill=trait_type)) +
    geom_violin(linewidth=0) +
    geom_boxplot(color="black", fill="white", linewidth=(0.5 / (ggplot2::.pt * 72.27/96)),
                 width=0.1, outlier.size = 0.5, outlier.shape = 16) +
    # geom_point(data = . %>% group_by(trait_type) %>% slice_max(n_signals), alpha=0.5, shape=16) +
    # geom_text_repel(data = . %>% group_by(trait_type) %>% slice_max(n_signals), aes(label=GWAS), size=2) +
    scale_fill_manual(values=trait_type_colors, guide="none") +
    scale_x_discrete(labels=c("cis-eQTLs" = "Cis-eQTLs",
                                "trans-eQTLs" = "Trans-eQTLs",
                                "all" = "Cis- and/or\ntrans-eQTLs")) +
    scale_y_continuous(name="GWAS signals colocalizing with eQTLs (%)") +
    facet_grid(cols=vars(trait_type), labeller=label_wrap_gen(width = 16, multi_line = TRUE)) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), axis.title.x = element_blank())

  p2

  ggsave(sprintf("colocalizations_cis_vs_trans_%s.pdf", Sys.Date()), p2, width=108, height=64, units="mm")


  colocalization_stats_exclusive <- joined %>%
    group_by(GWAS) %>%
    summarise(
      n_signals = n(),
      n_coloc_both = sum(cis & trans),
      n_coloc_cis = sum(cis & !trans),
      n_coloc_trans = sum(trans & !cis),
    ) %>%
    pivot_longer(starts_with("n_coloc"),
                 names_pattern="n_coloc_(.+)",
                 values_to="n_coloc", names_to="eqtl_type") %>%
    mutate(
      trait_type = trait_categories[GWAS],
      disease_trait = GWAS %in% disease_traits,
      percentage = n_coloc / n_signals * 100,
      eqtl_type = ordered(eqtl_type, labels=c("Only cis-eQTLs", "Only trans-eQTLs", "Both", "All"), levels=c("cis", "trans", "both", "all"))) %>%
    left_join(gwas_selection, by = c("GWAS" = "trait_identifier")) %>%
    select(trait_name, everything()) %>%
    group_by(trait_name) %>%
    mutate(summed_n_coloc = sum(n_coloc)) %>%
    ungroup() %>%
    arrange(desc(summed_n_coloc)) %>%
    mutate(trait_name = ordered(trait_name, levels=rev(unique(trait_name))))

  colocalization_stats_proc <- colocalization_stats_exclusive %>%
    group_by(GWAS, trait_name, trait_type, n_signals) %>%
    arrange(eqtl_type) %>%
    summarise(n_coloc = sum(n_coloc)) %>%
    mutate(n_no_coloc = n_signals - n_coloc) %>%
    pivot_longer(cols = c(n_coloc, n_no_coloc), names_to = "type", values_to = "n") %>%
    mutate(type = ordered(type, labels = c("Not colocalizing signals", "Colocalizing signals"), levels = c("n_no_coloc", "n_coloc")))

  colocalization_stats_exclusive_summarised <- colocalization_stats_exclusive %>%
    group_by(GWAS, trait_name, trait_type, n_signals, disease_trait) %>%
    arrange(desc(eqtl_type)) %>%
    summarise(
      n_both = sum(n_coloc[eqtl_type == "Both"]),
      n_total = sum(n_coloc),
      coloc_lab = paste(n_coloc, sep = " / ", collapse = " / "),
      prop_both = n_coloc[eqtl_type == "Both"] / sum(n_coloc),
      sum_coloc = sum(n_coloc))

  model_disease_traits <- glm(
    cbind(n_both, n_total - n_both) ~ disease_trait,
    family = binomial,
    data = colocalization_stats_exclusive_summarised
  )

  summary(model_disease_traits)

  model_trait_types <- glm(
    cbind(n_both, n_total - n_both) ~ trait_type,
    family = binomial,
    data = colocalization_stats_exclusive_summarised
  )

  summary(model_trait_types)

  # Trait type coefficients
  trait_type_coef <- tidy(model_trait_types) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      trait_type = gsub("trait_type", "", term),
      label = paste0(
        "OR=", round(exp(estimate), 2),
        "\np=", format.pval(p.value, digits=2)
      )
    )

  # Disease coefficients
  disease_coef <- tidy(model_disease_traits) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      disease_trait = as.logical(gsub("disease_trait", "", term)),
      label = paste0(
        "OR=", round(exp(estimate), 2),
        "\np=", format.pval(p.value, digits=2)
      )
    )

  trait_type_annotation <- colocalization_stats_exclusive_summarised %>%
    group_by(trait_type) %>%
    summarise(
      y=max(prop_both, na.rm=TRUE),
      .groups="drop"
    ) %>%
    left_join(trait_type_coef, by="trait_type") %>%
    mutate(
      y=y + 0.05
    )


  disease_annotation <- colocalization_stats_exclusive_summarised %>%
    group_by(disease_trait) %>%
    summarise(
      y=max(prop_both, na.rm=TRUE),
      .groups="drop"
    ) %>%
    left_join(disease_coef, by="disease_trait") %>%
    mutate(
      y=y + 0.05
    )


  trait_type_prop_plot <- ggplot(colocalization_stats_exclusive_summarised, aes(trait_type, prop_both)) +
    geom_violin() +
    geom_beeswarm(aes(size = sum_coloc, color=trait_type), shape=16, cex=3, alpha=0.5) +
    scale_color_manual(values=trait_type_colors) +
    geom_text(
      data=trait_type_annotation,
      aes(
        x=trait_type,
        y=y,
        label=label
      ),
      inherit.aes=FALSE,
      size=6/.pt
    ) +
    ggtitle(
      label="Difference in proportion for trait types\ncompared to blood-cell-type-composition traits",
    ) +
    xlab("Trait type") +
    ylab(
      "Proportion of eQTL colocalizing GWAS CSs colocalizing to\ncis+trans eQTLs compared to just cis-eQTLs or just trans-eQTLs"
    ) +
    coord_cartesian(clip="off") +
    theme(
      legend.position="null",
      plot.margin=margin(5,5,40,5)
    )

  disease_prop_plot <- ggplot(colocalization_stats_exclusive_summarised, aes(disease_trait, prop_both)) +
    geom_violin() +
    geom_beeswarm(aes(size = sum_coloc, color=trait_type), shape=16, cex=3, alpha=0.5) +
    scale_color_manual(values=trait_type_colors) +
    geom_text(
      data=disease_annotation,
      aes(
        x=disease_trait,
        y=y,
        label=label
      ),
      inherit.aes=FALSE,
      size=6/.pt
    ) +
    ggtitle(
      "Difference in proportion for disease\ntraits vs non-disease traits"
    ) +
    xlab("Disease trait") +
    ylab(
      "Proportion of eQTL colocalizing GWAS CSs colocalizing to\ncis+trans eQTLs compared to just cis-eQTLs or just trans-eQTLs"
    ) +
    coord_cartesian(clip="off") +
    theme(
      legend.position="right",
      plot.margin=margin(5,5,40,5)
    )

  combined_glm <- trait_type_prop_plot + disease_prop_plot + plot_layout(ncol=2, axes="collect", widths=c(7, 3))
  ggsave("glm_combined.pdf", combined_glm, width=180, height=120, units = "mm")

  p_n_signal_coloc <- ggplot(colocalization_stats_proc, aes(n, trait_name, fill = type)) +
    geom_col(position = position_stack()) +
    scale_fill_manual(values = c("Colocalizing signals" = "#D55E00", "Not colocalizing signals" = "grey70")) +
    xlab("Number of credible sets") +
    facet_wrap(vars(trait_type), scales="free_y", space="free_y", ncol=1) +
    theme(panel.grid.major.x = element_line(color="grey70"),
          strip.text = element_text(hjust = 0), legend.position = "bottom", axis.title.y = element_blank())

  p_prop_signal_coloc_type <- ggplot(colocalization_stats_exclusive, aes(n_coloc, trait_name, fill=eqtl_type)) +
    geom_col(position = position_fill()) +
    facet_wrap(vars(trait_type), scales="free_y", space="free_y", ncol=1) +
    xlab("Proportion of colocalizations") +
    theme(panel.grid.major.x = element_line(color="grey70"),
          strip.text = element_text(hjust = 0), legend.position = "bottom",
          axis.text.y = element_blank(), axis.title.y = element_blank())

  p_count_signal_coloc_type <- ggplot(colocalization_stats_exclusive_summarised, aes(0, trait_name, label = coloc_lab)) +
    geom_text(position = position_fill(), size = 6/.pt, hjust=0) +
    facet_wrap(vars(trait_type), scales="free_y", space="free_y", ncol=1) +
    xlab("Counts") +
    theme(strip.text = element_blank(), legend.position = "bottom",
          axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.line = element_blank(), axis.ticks = element_blank(),
    )

  combined <- p_n_signal_coloc + p_prop_signal_coloc_type + p_count_signal_coloc_type + plot_layout(widths = c(3, 3, 1))

  ggsave("combined_rev2.pdf", combined, width=180, height=240, units="mm")
  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}