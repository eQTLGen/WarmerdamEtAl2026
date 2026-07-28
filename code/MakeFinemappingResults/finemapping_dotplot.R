#!/usr/bin/env Rscript


# Load libraries
library(data.table)
library(tidyverse)
library(ggrepel)
library(arrow)
library(ggpubr)
library(IGUtilityPackage)

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

identify_lead_eqtls <- function(sig_effects, p_thresh=5e-8, window = 1e6) {

  lead_snps <- sig_effects %>%
    group_by(phenotype) %>%
    group_modify(~ IdentifyLeadSNPs(
      .x, Pthresh = as.numeric(p_thresh), window = as.numeric(window),
      snp_id_col = "variant", snp_chr_col = "chromosome", snp_pos_col = "bp", beta_col = "beta", se_col = "standard_error"))

  lead_snps0 <- sig_effects[paste0(sig_effects$variant, sig_effects$phenotype) %in% paste0(lead_snps$SNP, lead_snps$phenotype), ]

  lead_snps0$type <- "trans"

  lead_snps0$seqid_numeric <- case_when(
    lead_snps0$seqid == "X" ~ 23, lead_snps0$seqid == "Y" ~ 24, lead_snps0$seqid == "XY" ~ 25, lead_snps0$seqid == "MT" ~ 26,
    TRUE ~ as.numeric(as.character(lead_snps0$seqid)))

  lead_snps0[lead_snps0$chromosome == lead_snps0$seqid_numeric &
               lead_snps0$bp > lead_snps0$start - window &
               lead_snps0$bp < lead_snps0$end + window, ]$type <- "cis"

  # In case of two cis-eQTL SNPs, take only strongest
  lead_snps3 <- lead_snps0 %>% filter(type == "cis") %>% group_by(phenotype) %>% mutate(cis_effect_nr = row_number(abs(beta/standard_error)))

  lead_snps2 <- bind_rows(lead_snps0[lead_snps0$type == "trans", ], lead_snps3)

  lead_snps2 <- lead_snps2 %>% arrange(chromosome, bp)

  return(lead_snps2)
}

# Declare function definitions

GetGenome <- function(build) {
  if (build == "hg19") {
    chromosome <- data.frame(chr = c("1", "2", "3", "4", "5", "6", "7", "X", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "20", "Y", "19", "22", "21"), size = c(2.49251e+08, 2.43199e+08, 1.98022e+08, 1.91154e+08, 1.80915e+08, 1.71115e+08, 1.59139e+08, 1.55271e+08, 1.46364e+08, 1.41213e+08, 1.35535e+08, 1.35007e+08, 1.33852e+08, 1.1517e+08, 1.0735e+08, 1.02531e+08, 9.03548e+07, 8.11952e+07, 7.80772e+07, 6.30255e+07, 5.93736e+07, 5.9129e+07, 5.13046e+07, 4.81299e+07))
  } else if (build == "hg38") {
    chromosome <- data.frame(chr = c("1", "2", "3", "4", "5", "6", "7", "X", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "20", "Y", "19", "22", "21"), size = c(2.48956e+08, 2.42194e+08, 1.98296e+08, 1.90215e+08, 1.81538e+08, 1.70806e+08, 1.59346e+08, 1.56041e+08, 1.45139e+08, 1.38395e+08, 1.35087e+08, 1.33797e+08, 1.33275e+08, 1.14364e+08, 1.07044e+08, 1.01991e+08, 9.03383e+07, 8.32574e+07, 8.03733e+07, 6.44442e+07, 5.86176e+07, 5.72274e+07, 5.08185e+07, 4.671e+07))
  } else {
    stop("build has to be one of hg19 or hg38!")
  }
  chromosome$chr <- as.character(chromosome$chr)
  xy <- chromosome[chromosome$chr %in% c("X", "Y"),]
  chromosome <- chromosome[chromosome$chr %in% c(1:22),]
  chromosome$chr <- as.numeric(chromosome$chr)
  chromosome$abs_coord <- 0
  chromosome <- chromosome[order(chromosome$chr),]
  chromosome[1,]$abs_coord <- chromosome[1,]$size
  for (i in 2:nrow(chromosome)) {
    chromosome[i,]$abs_coord <- chromosome[i - 1,]$abs_coord + chromosome[i,]$size
  }
  xy$abs_coord <- NA
  chromosome <- rbind(chromosome, xy)
  chromosome$abs_coord[23] <- chromosome$abs_coord[22] + chromosome$size[23]
  chromosome$abs_coord[24] <- chromosome$abs_coord[23] + chromosome$size[24]
  chromosome$tick_mark <- chromosome$abs_coord - chromosome$size / 2
  chromosome$tick_mark[1] <- chromosome$abs_coord[1] / 2
  return(chromosome)
}


GetChromosomeColors <- function(chromosome, base, lighten_even = 0.1) {
  # Create alternating pattern (1 = odd, 2 = even)
  chromosome_pattern <- rep_len(c(1, 2), length.out = nrow(chromosome))
  names(chromosome_pattern) <- chromosome$chr

  result <- c()

  for (base_name in names(base)) {
    base_color <- base[[base_name]]
    for (chr_name in names(chromosome_pattern)) {
      chr_type <- chromosome_pattern[[chr_name]]
      color <- if (chr_type == 2) {
        lighten(base_color, lighten_even)
      } else {
        base_color
      }
      result_name <- paste(base_name, chr_name, sep="_chr")
      result[[result_name]] <- color
    }
  }

  return(result)
}

PlotDotPlot <- function(InputAssoc, genesWindow, pvaluecol, zscorecol = NULL, snpchrcol,
                        snpchrposcol, genechrcol, genechrposcol, type = "P", colorcol = NULL,
                        alpha = 0.15, sideplots = FALSE, window = 1e+06, build = "hg19", color_scale=NULL)
{
  chromosome <- GetGenome(build)

  p <- ggplot(chromosome, aes(x = abs_coord, y = abs_coord)) +
    theme_bw() +
    geom_vline(xintercept = c(chromosome$abs_coord[1:22],
                              0), colour = "lightgrey") +
    geom_hline(yintercept = c(chromosome$abs_coord,
                              0), colour = "lightgrey") +
    scale_x_continuous(breaks = chromosome$tick_mark,
                       labels = paste("chr", chromosome$chr, sep = ""),
                       limits = c(0, max(chromosome$abs_coord[1:22])),
                       expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = chromosome$tick_mark,
                       labels = paste("chr", chromosome$chr, sep = ""),
                       limits = c(0, max(chromosome$abs_coord)),
                       expand = c(0.01, 0.01)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab(paste0("variant position (",
                build, ")")) +
    ylab(paste0("gene center position (",
                build, ")"))

  trans <- InputAssoc
  if (ncol(trans) < 5) {
    stop("Input seems to contain too few columns!")
  }
  if (!pvaluecol %in% colnames(trans)) {
    stop("P-value column not found from input!")
  }
  if (!snpchrcol %in% colnames(trans)) {
    stop("SNP chr column not found from input!")
  }
  if (!snpchrposcol %in% colnames(trans)) {
    stop("SNP chr position column not found from input!")
  }
  if (!genechrcol %in% colnames(trans)) {
    stop("Gene chr column not found from input!")
  }
  if (!genechrposcol %in% colnames(trans)) {
    stop("Gene chr position column not found from input!")
  }
  if (!is.null(colorcol)) {
    if (!colorcol %in% colnames(trans)) {
      stop("Color column is specified but not present in data!")
    }
  }
  trans$abs_SNPChrPos <- 0
  trans[trans[[snpchrcol]] == 1,]$abs_SNPChrPos <- trans[trans[[snpchrcol]] ==
                                                           1,][[snpchrposcol]]
  for (i in 2:22) {
    trans[trans[[snpchrcol]] == i,]$abs_SNPChrPos <- trans[trans[[snpchrcol]] ==
                                                             i,][[snpchrposcol]] + chromosome[chromosome$chr ==
                                                                                                i - 1,]$abs_coord
  }
  trans$abs_ProbeChrPos <- 0
  trans[trans[[genechrcol]] == 1,]$abs_ProbeChrPos <- trans[trans[[genechrcol]] ==
                                                              1,][[genechrposcol]]
  trans$abs_ProbeChrPos <- as.numeric(trans$abs_ProbeChrPos)
  for (i in c(2:24)) {
    trans[trans[[genechrcol]] == chromosome$chr[i],]$abs_ProbeChrPos <- (
      as.numeric(trans[trans[[genechrcol]] == chromosome$chr[i],][[genechrposcol]]) +
        as.numeric(chromosome[chromosome$chr == chromosome$chr[i - 1],]$abs_coord))
  }
  trans$Pvalue <- trans[[pvaluecol]]
  if (!is.null(zscorecol)) {
    trans$OverallZScore <- trans[[zscorecol]]
  }
  message(paste(nrow(trans), "rows in QTL file."))

  RollingSum <- function(x, window = 1e+06) {
    Sum <- rep(0, times = length(x))
    for (i in 1:length(x)) {
      Sum[i] <- length(x[x > x[i] - window & x < x[i] +
        window])
    }
    return(Sum)
  }

  trans <- trans[order(trans$abs_SNPChrPos),]
  trans$SnpSum <- RollingSum(trans$abs_SNPChrPos, window = window)
  trans <- trans[order(trans$abs_ProbeChrPos),]
  trans$GeneSum <- RollingSum(trans$abs_ProbeChrPos, window = window)
  trans$SnpChr <- trans[[snpchrcol]]
  trans$GeneChr <- trans[[genechrcol]]
  trans$SnpChr <- factor(trans$SnpChr, levels = c(1:22,
                                                  "X", "Y"))
  trans$GeneChr <- factor(trans$GeneChr, levels = c(1:22,
                                                    "X", "Y"))

  if (!is.null(colorcol)) {
    trans$temp <- trans[[colorcol]]
    trans$temp_chr_col <- paste(trans$temp, trans$GeneChr, sep="_chr")
  }

  if (type == "P") {
    if (!is.null(colorcol)) {

      point_plot <- p +
        geom_point(data = as.data.frame(trans),
                   aes(x = abs_SNPChrPos, y = abs_ProbeChrPos, size = -log10(Pvalue), colour = temp),
                   alpha = alpha, shape=16) +
        scale_size_continuous(breaks = c(10, 50, 100,
                                         200, 300), range = c(1, 5), guide = guide_legend(title = expression(paste(-log[10]("P-value"))))) +
        ylab(paste0("Gene position (", build, ")")) +
        scale_color_manual(values=color_scale, guide=NULL)
    }
    else {
      point_plot <- p +
        geom_point(data = as.data.frame(trans),
                   aes(x = abs_SNPChrPos, y = abs_ProbeChrPos, size = -log10(Pvalue)),
                   alpha = alpha, colour = "black", shape=16) +
        scale_size_continuous(breaks = c(10,
                                         50, 100, 200, 300), range = c(1, 5), guide = guide_legend(title = expression(paste(-log[10]("P-value"))))) +
        ylab(paste0("Gene position (", build, ")"))
    }
  }
  if (type == "Z") {
    if (!is.null(colorcol)) {
      point_plot <- p +
        geom_point(data = as.data.frame(trans),
                   aes(x = abs_SNPChrPos, y = abs_ProbeChrPos, size = -log10(OverallZScore), colour = temp),
                   alpha = alpha, shape=16) +
        scale_size_continuous(breaks = c(10, 50, 100,
                                         200, 300), range = c(1, 5), guide = guide_legend(title = expression(paste(-log[10]("P-value"))))) +
        ylab(paste0("Gene position (", build, ")")) +
        scale_color_manual(values=color_scale)
    }
    else {
      point_plot <- p +
        geom_point(data = as.data.frame(trans),
                   aes(x = abs_SNPChrPos, y = abs_ProbeChrPos, size = abs(OverallZScore)),
                   alpha = alpha, colour = "black", shape=16) +
        scale_size_continuous(range = c(1,
                                        5), guide = guide_legend(title = expression(paste("meta-analysis Z-score")))) +
        ylab(paste0("Gene position (", build, ")"))
    }
  }
  if (sideplots == TRUE) {

    GenewiseSummary <- trans %>%
      group_by(phenotype, gene_name, abs_ProbeChrPos, GeneChr) %>%
      summarise(n_signals = n()) %>%
      mutate(label = sprintf("%s (%s)", ifelse(is.na(gene_name), phenotype, gene_name), n_signals)) %>%
      ungroup()

    # Use
    locusInputDt <- as.data.table(trans[,c("variant", "chromosome", "bp", "SnpChr", "abs_SNPChrPos")])
    locusInputDt$bp_end <- locusInputDt$bp

    # Set keys
    setkey(genesWindow, seqid, window_start, window_end)
    setkey(locusInputDt, chromosome, bp, bp_end)

    # Perform overlap join
    overlaps <- foverlaps(locusInputDt, genesWindow, nomatch = 0)

    LocuswiseSummary <- overlaps %>%
      group_by(gene_id, gene_name, window_start, window_end, SnpChr) %>%
      summarise(
        n_variants=n(),
        variant_list = list(variant),
        abs_SNPChrPos = (min(abs_SNPChrPos)+max(abs_SNPChrPos))/2) %>%
      mutate(
        n_variants_capped = pmin(500, n_variants),
        capped = n_variants!=n_variants_capped,
        label = sprintf("%s (%s)", ifelse(is.na(gene_name), gene_id, gene_name), n_variants)) %>%
      ungroup()

    p1 <- ggplot(LocuswiseSummary, aes(y = n_variants_capped, x = abs_SNPChrPos,
                                       colour = SnpChr)) +
      theme_bw() +
      geom_point(aes(shape=capped)) +
      geom_segment(
        data = . %>% filter(n_variants > 500), aes(x = abs_SNPChrPos, y = n_variants_capped, yend = max(LocuswiseSummary$n_variants_capped) * 1.2)) +
      geom_text_repel(
        data = . %>% filter(n_variants > 500), aes(abs_SNPChrPos, max(LocuswiseSummary$n_variants_capped) * 1.2, label = label),
        nudge_y=1, hjust = "left", direction="x", min.segment.length = 0, angle=90, color="black", size=3.88*0.8) +
      geom_vline(xintercept = c(chromosome$abs_coord[1:22],
                                0), colour = "lightgrey", alpha = 0.3) +
      scale_x_continuous(breaks = chromosome$tick_mark,
                         labels = paste("chr", chromosome$chr, sep = ""),
                         limits = c(0, max(chromosome$abs_coord[1:22])),
                         expand = c(0.01, 0.01)) +
      theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title.x = element_blank()) +
      scale_shape_manual(values = c(16, 17)) +
      scale_colour_manual(values = c(`1` = "darkgrey",
                                     `3` = "darkgrey", `5` = "darkgrey", `7` = "darkgrey",
                                     `9` = "darkgrey", `11` = "darkgrey", `13` = "darkgrey",
                                     `15` = "darkgrey", `17` = "darkgrey", `19` = "darkgrey",
                                     `21` = "darkgrey", X = "darkgrey", `2` = "black",
                                     `4` = "black", `6` = "black", `8` = "black",
                                     `10` = "black", `12` = "black", `14` = "black",
                                     `16` = "black", `18` = "black", `20` = "black",
                                     `22` = "black", Y = "black")) +
      xlab(paste0("SNP position (",
                  build, ")")) +
      scale_y_continuous(
        breaks = seq(from = 0,
                     to = 500,
                     by = 200),
        limits = c(0, 500 * 3),
        name="N primary eQTLs\n±1Mb of gene") +
      guides(colour = "none")

    p2 <- ggplot(GenewiseSummary, aes(y = n_signals, x = abs_ProbeChrPos,
                                      colour = GeneChr)) +
      theme_bw() +
      geom_vline(xintercept = c(chromosome$abs_coord,
                                0), colour = "lightgrey", alpha = 0.3) +
      geom_point() +
      geom_segment(data = . %>% slice_max(n_signals, n=10), aes(x = abs_ProbeChrPos, y = n_signals, yend = max(GenewiseSummary$n_signals) * 1.2)) +
      geom_text_repel(
        data = . %>% slice_max(n_signals, n=10), aes(abs_ProbeChrPos, max(GenewiseSummary$n_signals) * 1.2, label = label),
        nudge_y=2, hjust = "left", direction="y", min.segment.length = 0, color="black", size = 3.88*0.8) +
      scale_x_continuous(breaks = chromosome$tick_mark,
                         labels = paste("chr", chromosome$chr, sep = ""),
                         limits = c(0, max(chromosome$abs_coord)), expand = c(0.01, 0.01)) +
      scale_colour_manual(values = c(`1` = "darkgrey",
                                     `3` = "darkgrey", `5` = "darkgrey", `7` = "darkgrey",
                                     `9` = "darkgrey", `11` = "darkgrey", `13` = "darkgrey",
                                     `15` = "darkgrey", `17` = "darkgrey", `19` = "darkgrey",
                                     `21` = "darkgrey", X = "darkgrey", `2` = "black",
                                     `4` = "black", `6` = "black", `8` = "black", `10` = "black",
                                     `12` = "black", `14` = "black", `16` = "black", `18` = "black",
                                     `20` = "black", `22` = "black", Y = "black")) +
      theme(axis.text.y = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      xlab(paste0("Gene position (", build, ")")) +
      coord_flip() +
      scale_y_continuous(
        breaks = seq(from = 0,
                     to = round(max(GenewiseSummary$n_signals), digits=-1),
                     by = round(max(GenewiseSummary$n_signals) / 3, digits=-1)),
        limits = c(0, max(round(GenewiseSummary$n_signals) * 3)),
        name = "N primary eQTLs per gene") +
      guides(colour = "none")

    layout <- "\n    E##\n    D##\n    ABC\n    "

    dotplot <- point_plot + plot_spacer() + p2 + plot_spacer() + p1 +
      plot_layout(design = layout, widths = c(5, 0.1, 2), heights = c(2, 0.1, 5), guides='collect') &
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position="bottom")

    return(dotplot)
  }
  else {
    dotplot <- point_plot
  }
  return(dotplot)
}


# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Define variables
  outcome_variables <- list()

  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))

  available_genes <- fread(variables["available_genes"], header=F)$V1
  all_significant <- fread(variables["all_significant"])

  # Process input
  variant_reference <- read_parquet(variables["variant_reference"])
  gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
    filter(type == "gene") %>% select(gene_id, gene_name, seqid, start, end, strand) %>%
    mutate(seqid_numeric = case_when(
      seqid == "X" ~ 23, seqid == "Y" ~ 24, seqid == "XY" ~ 25, seqid == "MT" ~ 26,
      TRUE ~ as.numeric(as.character(seqid))))

  outcome_variables["n_available_genes_meta_analysis_in_gene_reference"] <- sum(available_genes %in% gene_reference$gene_id)

  lead_snps2 <- identify_lead_eqtls(all_significant)

  fwrite(lead_snps2, "../../lead_variants.txt.gz", sep="\t", col.names=T, row.names=F)

  significant_lead_snps <- deframe(lead_snps2 %>%
                                     summarise(
                                       n_genes_assessed = sum(available_genes %in% gene_reference$gene_id),
                                       n_egenes = n_distinct(phenotype),
                                       n_eqtls = n(),
                                       n_eqtls_cis = sum(type=="cis"),
                                       n_eqtls_trans = sum(type=="trans"),
                                       n_egenes_cis = n_distinct(phenotype[type=="cis"]),
                                       n_egenes_trans = n_distinct(phenotype[type=="trans"]),
                                       n_variants = n_distinct(variant),
                                       p_egenes = n_egenes/n_genes_assessed*100,
                                       p_egenes_cis = n_egenes_cis/n_genes_assessed*100,
                                       p_egenes_trans = n_egenes_trans/n_genes_assessed*100) %>% pivot_longer(everything()))

  outcome_variables[["signficant_eqtls_per_2Mb"]] <- significant_lead_snps

  # Window sizes
  cis_window <- 1e6
  twilight_window <- 5e6

  outcome_variables["cis_window"] <- cis_window
  outcome_variables["twilight_window"] <- twilight_window

  # Windows to plot
  window_size <- 1e6
  window <- c("start"=-window_size/2, "end"=window_size/2, "by" = 2e5)

  # Color pallete
  gene_locus_colors <- c("upstream" = "#009E73", "intragenic" = "#D55E00", "downstream" = "#0072B2")
  cis_trans_color_scale <- c("cis" = "#D7817F", "trans" = "#0072B2")

  finemapped_eqtls <- fread(variables["independent_variant_set"]) %>%
    mutate(Z = beta / standard_error,
           absZ = abs(Z),
           absBeta = abs(beta),
           p_value = 2*pnorm(absZ, lower.tail=F),
           variant_class = case_when(
             str_detect(variant, "HGSV") ~ "SV",
             str_length(eff_allele) > str_length(non_eff_allele) ~ "Insertion",
             str_length(eff_allele) < str_length(non_eff_allele) ~ "Deletion",
             str_length(eff_allele) == str_length(non_eff_allele) ~ "SNP"))

  build <- "hg38"

  dotplot <- PlotDotPlot(as_tibble(finemapped_eqtls),
                         genesWindow = NA,
                         pvaluecol = "p_value",
                         snpchrcol = "chromosome",
                         snpchrposcol = "bp",
                         genechrcol = "seqid",
                         genechrposcol = "start",
                         sideplots = F,
                         build = build,
                         zscorecol = "Z",
                         colorcol= "type", color_scale=cis_trans_color_scale)

  dotplot

  ggsave("dotplot.pdf", dotplot, width=180, height=144, units = "mm")

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}