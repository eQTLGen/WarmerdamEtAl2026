library(data.table)
library(IGUtilityPackage)
library(tidyverse)
library(qvalue)
library(patchwork)
library(ggrepel)
library(scales)
library(progressr)
library(tictoc)
library(ggrepel)
library(grid)
library(gridGraphics)

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

# Rb calculations were done in separate Julia script, hereby the visualisation

variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))

rb1 <- fread("rb_results/cis_trans_replication_rb_ClumpedVariants_2025-12-17.tsv")
rb2 <- fread("rb_results/ccTrans_intracellularTrans_replication_rb_ClumpedVariants_2025-12-17.tsv")
eqtl_counts <- fread("rb_results/eQTL_counts_rb_ClumpedVariants_2025-12-17.tsv")
rb <- rbind(rb1, rb2)


rb[type == "cis"]$type <- "cis-eQTLs"
rb[type == "trans"]$type <- "trans-eQTLs"
rb[type == "ct"]$type <- "cell-type-composition trans-eQTLs"
rb[type == "intracellular"]$type <- "intracellular trans-eQTLs"

# Add Rb counts and fractions, how large fraction was tested in each replication dataset

colnames(eqtl_counts)[1] <- "dataset"
eqtl_counts[type == "cis"]$type <- "cis-eQTLs"
eqtl_counts[type == "trans"]$type <- "trans-eQTLs"
eqtl_counts[type == "ct"]$type <- "cell-type-composition trans-eQTLs"
eqtl_counts[type == "intracellular"]$type <- "intracellular trans-eQTLs"

finemapped <- fread("rb_results/replications_formatted_for_rb_ClumpedVariants_2025-12-17.txt")
finemapped <- finemapped[, c(1, 2), with = FALSE]
finemapped$trans_type <- NA
ct <- fread(variables["blood_gwas_coloc"])
ct <- paste(ct$variant, ct$phenotype, sep = "-")
finemapped <- finemapped[type != "intermediate"]
finemapped <- as.data.frame(finemapped)

finemapped[!finemapped$eQTL %in% ct & !finemapped$type %in% "cis", ]$trans_type <- "intracellular trans-eQTL"
finemapped[finemapped$eQTL %in% ct & !finemapped$type %in% "cis", ]$trans_type <- "cell-type-composition trans-eQTL"

finemapped <- as.data.table(finemapped)
nr_cis_eQTLs <- length(unique(finemapped[type == "cis"]$eQTL))
nr_trans_eQTLs <- length(unique(finemapped[type == "trans"]$eQTL))
nr_cc_eQTLs <- length(unique(finemapped[trans_type == "cell-type-composition trans-eQTL"]$eQTL))
nr_intracellular_eQTLs <- length(unique(finemapped[trans_type == "intracellular trans-eQTL"]$eQTL))

eqtl_counts$discovery_eQTLs <- NA
eqtl_counts <- as.data.frame(eqtl_counts)
eqtl_counts[eqtl_counts$type %in% "cis-eQTLs", ]$discovery_eQTLs <- nr_cis_eQTLs
eqtl_counts[eqtl_counts$type %in% "trans-eQTLs", ]$discovery_eQTLs <- nr_trans_eQTLs
eqtl_counts[eqtl_counts$type %in% "cell-type-composition trans-eQTLs", ]$discovery_eQTLs <- nr_cc_eQTLs
eqtl_counts[eqtl_counts$type %in% "intracellular trans-eQTLs", ]$discovery_eQTLs <- nr_intracellular_eQTLs

eqtl_counts <- as.data.table(eqtl_counts)
eqtl_counts$fraction_eQTLs_tested <- eqtl_counts$nr_eQTLs / eqtl_counts$discovery_eQTLs

rb <- merge(rb, eqtl_counts, by = c("dataset", "type"))

rb$significance <- "≥0.05"
rb[p<0.05]$significance <- "<0.05"

rb$dataset_type <- "GTEx tissues bulk"
rb[str_detect(dataset, "scMeta")]$dataset_type <- "single-cell blood"
rb[str_detect(dataset, "scBrain")]$dataset_type <- "single-cell brain"
rb[str_detect(dataset, "pqtl")]$dataset_type <- "plasma pQTL"
rb[str_detect(dataset, "metaBrain")]$dataset_type <- "brain bulk"
rb[str_detect(dataset, "LCL")]$dataset_type <- "LCL bulk"

rb[type == "ct"]$type <- "cell-type-composition trans-eQTLs"
rb[type == "intracellular"]$type <- "intracellular trans-eQTLs"
rb[type == "cis"]$type <- "cis-eQTLs"
rb[type == "trans"]$type <- "trans-eQTLs"

rb$type <- factor(rb$type, levels = c("cis-eQTLs", "trans-eQTLs", "cell-type-composition trans-eQTLs", "intracellular trans-eQTLs"))

# Clean tissues
rb$dataset <- str_replace_all(rb$dataset, "_", " ")
rb$dataset <- str_replace(rb$dataset, "gtex ", "GTEx: ")
rb$dataset <- str_replace(rb$dataset, "scMeta ", "blood sc-eQTL: ")
rb$dataset <- str_replace(rb$dataset, "scBrain ", "brain sc-eQTL: ")
rb$dataset <- str_replace(rb$dataset, "pqtl", "pQTL")
rb$dataset <- str_replace(rb$dataset, "metaBrain", "brain bulk eQTL")
rb$dataset <- str_replace(rb$dataset, "meta LCL", "LCL bulk eQTL")

rb$dataset <- factor(rb$dataset, levels = rev(unique(rb$dataset)))

rb$indik <- "trans"
rb[type == "cis-eQTLs"]$indik <- "cis"

fwrite(rb, "rb_results/rb_results_cleaned_2025-12-17.tsv", sep = "\t")

p1 <- ggplot(rb[dataset_type == "GTEx tissues bulk"], aes(y = dataset, x = rb, color = type, alpha = significance)) +
  geom_vline(xintercept = c(seq(from = -1, to = -0.2, by = 0.2), seq(from = 0.2, to = 1, by = 0.2)), colour = "grey87", linetype = 1, size = 0.2) +
  geom_vline(xintercept = 0, colour = "black", linetype = 1) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(xmin = rb - rb_se, xmax = rb + rb_se),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  labs(y = "Dataset", x = "Rb", color = "Type") +
  scale_x_continuous(breaks = seq(-1, 1, 0.2)) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_colour_brewer(palette = "Set1") +
  scale_alpha_manual(values = c("≥0.05" = 0.3, "<0.05" = 1)) +
  theme(         # remove legend
    axis.text.x = element_blank(),     # remove x-axis text
    axis.ticks.x = element_blank(),
    legend.position="none"# remove x-axis ticks
  ) +
  labs(x = NULL, y = NULL) +
  facet_grid(~ indik)

p2 <- ggplot(rb[dataset_type == "LCL bulk"], aes(y = dataset, x = rb, color = type, alpha = significance)) +
  geom_vline(xintercept = c(seq(from = -1, to = -0.2, by = 0.2), seq(from = 0.2, to = 1, by = 0.2)), colour = "grey87", linetype = 1, size = 0.2) +
  geom_vline(xintercept = 0, colour = "black", linetype = 1) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(xmin = rb - rb_se, xmax = rb + rb_se),
                 width = 0.2,
                 position = position_dodge(width = 0.5)) +
  labs(y = "Dataset", x = "Rb", color = "Type") +
  scale_x_continuous(breaks = seq(-1, 1, 0.2)) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_colour_brewer(palette = "Set1") +
  scale_alpha_manual(values = c("≥0.05" = 0.3, "<0.05" = 1)) +
  theme(         # remove legend
    axis.text.x = element_blank(),     # remove x-axis text
    axis.ticks.x = element_blank(),
    legend.position="none"# remove x-axis ticks
  ) +
  labs(x = NULL, y = NULL) +
  facet_grid(~ indik)

p3 <- ggplot(rb[dataset_type == "single-cell blood"], aes(y = dataset, x = rb, color = type, alpha = significance)) +
  geom_vline(xintercept = c(seq(from = -1, to = -0.2, by = 0.2), seq(from = 0.2, to = 1, by = 0.2)), colour = "grey87", linetype = 1, size = 0.2) +
  geom_vline(xintercept = 0, colour = "black", linetype = 1) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(xmin = rb - rb_se, xmax = rb + rb_se),
                 width = 0.2,
                 position = position_dodge(width = 0.5)) +
  labs(y = "Dataset", x = "Rb", color = "Type") +
  scale_x_continuous(breaks = seq(-1, 1, 0.2)) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_colour_brewer(palette = "Set1") +
  scale_alpha_manual(values = c("≥0.05" = 0.3, "<0.05" = 1)) +
  theme(         # remove legend
    axis.text.x = element_blank(),     # remove x-axis text
    axis.ticks.x = element_blank(),
    legend.position="none"# remove x-axis ticks
  ) +
  labs(x = NULL, y = NULL) +
  facet_grid(~ indik)

p4 <- ggplot(rb[dataset_type == "brain bulk"], aes(y = dataset, x = rb, color = type, alpha = significance)) +
  geom_vline(xintercept = c(seq(from = -1, to = -0.2, by = 0.2), seq(from = 0.2, to = 1, by = 0.2)), colour = "grey87", linetype = 1, size = 0.2) +
  geom_vline(xintercept = 0, colour = "black", linetype = 1) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(xmin = rb - rb_se, xmax = rb + rb_se),
                 width = 0.2,
                 position = position_dodge(width = 0.5)) +
  labs(y = "Dataset", x = "Rb", color = "Type") +
  scale_x_continuous(breaks = seq(-1, 1, 0.2)) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_colour_brewer(palette = "Set1") +
  scale_alpha_manual(values = c("≥0.05" = 0.3, "<0.05" = 1)) +
  theme(         # remove legend
    axis.text.x = element_blank(),     # remove x-axis text
    axis.ticks.x = element_blank(),
    legend.position="none"# remove x-axis ticks
  ) +
  labs(x = NULL, y = NULL) +
  facet_grid(~ indik)

p5 <- ggplot(rb[dataset_type == "single-cell brain"], aes(y = dataset, x = rb, color = type, alpha = significance)) +
  geom_vline(xintercept = c(seq(from = -1, to = -0.2, by = 0.2), seq(from = 0.2, to = 1, by = 0.2)), colour = "grey87", linetype = 1, size = 0.2) +
  geom_vline(xintercept = 0, colour = "black", linetype = 1) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(xmin = rb - rb_se, xmax = rb + rb_se),
                 width = 0.2,
                 position = position_dodge(width = 0.5)) +
  labs(y = "Dataset", x = "Rb", color = "Type") +
  scale_x_continuous(breaks = seq(-1, 1, 0.2)) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_colour_brewer(palette = "Set1") +
  scale_alpha_manual(values = c("≥0.05" = 0.3, "<0.05" = 1)) +
  theme(         # remove legend
    axis.text.x = element_blank(),     # remove x-axis text
    axis.ticks.x = element_blank(),
    legend.position="none"# remove x-axis ticks
  ) +
  labs(x = NULL, y = NULL) +
  facet_grid(~ indik)

p6 <- ggplot(rb[dataset_type == "plasma pQTL"], aes(y = dataset, x = rb, color = type, alpha = significance)) +
  geom_vline(xintercept = c(seq(from = -1, to = -0.2, by = 0.2), seq(from = 0.2, to = 1, by = 0.2)), colour = "grey87", linetype = 1, size = 0.2) +
  geom_vline(xintercept = 0, colour = "black", linetype = 1) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(xmin = rb - rb_se, xmax = rb + rb_se),
                 width = 0.2,
                 position = position_dodge(width = 0.5)) +
  labs(y = "Dataset", x = "Rb", color = "Type") +
  scale_x_continuous(breaks = seq(-1, 1, 0.2)) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_colour_brewer(palette = "Set1") +
  scale_alpha_manual(values = c("≥0.05" = 0.3, "<0.05" = 1)) +
  labs(x = NULL, y = NULL) +
  theme(legend.position="none") +
  facet_grid(~ indik)

# Counts

rb2 <- melt(rb, id.vars = c("dataset", "type", "rb", "rb_se", "p", "significance", "dataset_type", "indik"))
rb2 <- rb2[variable %in% c("nr_eQTLs", "nr_genes", "nr_variants")]

p1_1 <- ggplot(rb[dataset_type == "GTEx tissues bulk"], aes(y = dataset, x = nr_genes, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_genes))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank()) +
  ggtitle("Genes")# remove x-axis ticks)
p1_2 <- ggplot(rb[dataset_type == "GTEx tissues bulk"], aes(y = dataset, x = nr_variants, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_variants))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank()) +
  ggtitle("Variants")# remove x-axis ticks)# remove x-axis ticks)
p1_3 <- ggplot(rb[dataset_type == "GTEx tissues bulk"], aes(y = dataset, x = nr_eQTLs, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_eQTLs))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank()) +
  ggtitle("eQTLs")# remove x-axis ticks)# remove x-axis ticks)


p2_1 <- ggplot(rb[dataset_type == "LCL bulk"], aes(y = dataset, x = nr_genes, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_genes))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank())# remove x-axis ticks)
p2_2 <- ggplot(rb[dataset_type == "LCL bulk"], aes(y = dataset, x = nr_variants, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_variants))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank())# remove x-axis ticks)
p2_3 <- ggplot(rb[dataset_type == "LCL bulk"], aes(y = dataset, x = nr_eQTLs, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_eQTLs))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank())# remove x-axis ticks)

p3_1 <- ggplot(rb[dataset_type == "single-cell blood"], aes(y = dataset, x = nr_genes, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_genes))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank())# remove x-axis ticks)
p3_2 <- ggplot(rb[dataset_type == "single-cell blood"], aes(y = dataset, x = nr_variants, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_variants))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank())# remove x-axis ticks)
p3_3 <- ggplot(rb[dataset_type == "single-cell blood"], aes(y = dataset, x = nr_eQTLs, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_eQTLs))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank())# remove x-axis ticks)

p4_1 <- ggplot(rb[dataset_type == "brain bulk"], aes(y = dataset, x = nr_genes, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_genes))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank())# remove x-axis ticks)
p4_2 <- ggplot(rb[dataset_type == "brain bulk"], aes(y = dataset, x = nr_variants, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_variants))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank())# remove x-axis ticks)
p4_3 <- ggplot(rb[dataset_type == "brain bulk"], aes(y = dataset, x = nr_eQTLs, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_eQTLs))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank())# remove x-axis ticks)

p5_1 <- ggplot(rb[dataset_type == "single-cell brain"], aes(y = dataset, x = nr_genes, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_genes))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank())# remove x-axis ticks)
p5_2 <- ggplot(rb[dataset_type == "single-cell brain"], aes(y = dataset, x = nr_variants, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_variants))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank())# remove x-axis ticks)
p5_3 <- ggplot(rb[dataset_type == "single-cell brain"], aes(y = dataset, x = nr_eQTLs, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_eQTLs))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),     # remove x-axis text
        axis.ticks.x = element_blank())# remove x-axis ticks)

p6_1 <- ggplot(rb[dataset_type == "plasma pQTL"], aes(y = dataset, x = nr_genes, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_genes))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position="none")# remove x-axis ticks)
p6_2 <- ggplot(rb[dataset_type == "plasma pQTL"], aes(y = dataset, x = nr_variants, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_variants))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position="none")# remove x-axis ticks)
p6_3 <- ggplot(rb[dataset_type == "plasma pQTL"], aes(y = dataset, x = nr_eQTLs, fill = type)) +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max(rb$nr_eQTLs))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(),     # remove x-axis text
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position="none")# remove x-axis ticks)

row1 <- (p1 | p1_1 | p1_2 | p1_3) + plot_layout(widths = c(10, 1, 1, 1))
row2 <- (p2 | p2_1 | p2_2 | p2_3) + plot_layout(widths = c(10, 1, 1, 1))
row3 <- (p3 | p3_1 | p3_2 | p3_3) + plot_layout(widths = c(10, 1, 1, 1))
row4 <- (p4 | p4_1 | p4_2 | p4_3) + plot_layout(widths = c(10, 1, 1, 1))
row5 <- (p5 | p5_1 | p5_2 | p5_3) + plot_layout(widths = c(10, 1, 1, 1))
row6 <- (p6 | p6_1 | p6_2 | p6_3) + plot_layout(widths = c(10, 1, 1, 1))

final_plot <- row1 / row2 / row3 / row4 / row5 / row6 +
  plot_layout(heights = c(7, 0.5, 0.8, 0.5, 1.25, 0.5))


ggsave("rb_results/rb_overview_ClumpedVariants.png",
       height = 10 * 1.2, width = 8 * 1, units = "in", dpi = 400)

ggsave("rb_results/rb_overview_ClumpedVariants.pdf",
       height = 10 * 1.2, width = 8 * 1, units = "in", dpi = 400)

ggsave("rb_results/rb_overview_ClumpedVariants.svg",
       height = 10 * 1.2, width = 8 * 1, units = "in", dpi = 400)
