#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)


# Declare constants


# Define functions
get_ld_block <- function(vars, ld_list) {

  # find which chromosomes contain these variants
  all_blocks <- list()

  for (chr in names(ld_list)) {

    R <- ld_list[[chr]]

    keep <- intersect(vars, rownames(R))

    if (length(keep) > 0) {
      R_sub <- R[keep, keep, drop = FALSE]
      all_blocks[[chr]] <- R_sub
    }
  }

  # block diagonal combine across chromosomes
  if (length(all_blocks) == 1) return(all_blocks[[1]])

  # combine safely
  bigR <- Matrix::bdiag(all_blocks)
  bigR <- as.matrix(bigR)

  # ensure row/col names
  rownames(bigR) <- colnames(bigR) <- unlist(lapply(all_blocks, rownames))

  bigR
}

compute_h2 <- function(vars, beta, chr, ld_list) {

  names(beta) <- vars

  if (length(beta) > 1) {
    # align order
    R <- ld_list[[chr]][vars, vars, drop = FALSE]

    # compute conditional beta's
    beta_joint_std_hat <- solve(R, beta)

    # Calculate joint heritability of all included variants
    return(as.numeric(t(beta_joint_std_hat) %*% R %*% beta_joint_std_hat))
  } else if (length(beta) == 1) {
    return(beta^2)
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
  # Define variables
  variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))
  location <- str_split(variables["available_genes"], "freeze3")[[1]][1]

  available_genes <- fread(variables["available_genes"], header=F)$V1
  gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
    filter(type == "gene") %>% select(gene_id, gene_name, seqid, start, end, strand) %>%
    mutate(TSS = if_else(strand == "+", start, end))

  # Get background set of genes
  # I.e. the maximum set of genes for which we could hypothetically detect eQTLs
  # - The intersect of the genes in the gene-reference file
  # - and in the set of ensembl IDs for which the meta-analysis was run
  gene_reference_filtered <- gene_reference %>% filter(gene_id %in% available_genes)
  tested_genes <- gene_reference_filtered %>% pull(gene_id)

  # Process input
  ldsc <- bind_rows(
    mapply(fread, Sys.glob(variables["ldsc_output"]), SIMPLIFY=F)) %>%
    mutate(h2_p = pnorm(abs(h2_obs/h2_obs_se), lower.tail = F) * 2,
           type = "polygenic") %>%
    rename(c("phenotype" = "gene_id")) %>%
    filter(phenotype %in% tested_genes)

  eur_genes <- ldsc %>% pull(phenotype)

  eqtls <- fread(variables["independent_variant_set"]) %>%
    select(variant_index, chromosome, phenotype, type, eaf, beta, eaf_eur, beta_eur) %>%
    filter(phenotype %in% eur_genes) %>%
    mutate(
      typeCisTrans = if_else(type == "cis", "cis", "trans"),
      beta_std = beta * sqrt(2 * eaf * (1 - eaf)),
      beta_eur_std = beta_eur * sqrt(2 * eaf_eur * (1 - eaf_eur)))
  # Process output

  groups <- eqtls %>% group_by(phenotype, typeCisTrans, chromosome) %>%
    group_split()

  ld_matrices <- readRDS(variables["finemapped_variant_ld"])

  results <- lapply(groups, function(df) {
    data.frame(
      phenotype = unique(df$phenotype),
      type = unique(df$typeCisTrans),
      chromosome = unique(df$chromosome),
      n_vars = nrow(df),
      h2_obs = compute_h2(as.character(df$variant_index), df$beta_eur_std, df$chromosome[1], ld_matrices)
    )
  })

  results_df <- bind_rows(results) %>% as_tibble()

  results_df2 <- results_df %>% group_by(phenotype, type) %>%
    summarise(n_vars = sum(n_vars), h2_obs = sum(h2_obs))

  h2_estimates <- bind_rows(ldsc, results_df2) %>%
    mutate(h2 = pmax(0, pmin(1, h2_obs))) %>%
    mutate(type = factor(type, levels=c("polygenic", "cis", "trans"), labels=c("polygenic", "cis-eQTLs", "trans-eQTLs"), ordered=T))

  h2_estimates_wide <- h2_estimates %>% select(h2_obs, h2_obs_se, phenotype, type) %>%
    pivot_wider(id_cols = "phenotype", values_from = c("h2_obs", "h2_obs_se"), names_from="type", names_sep="_", values_fill=0) %>%
    mutate(`eQTL total h2` = `h2_obs_cis-eQTLs` + `h2_obs_trans-eQTLs`, `Total common variant heritability` = `eQTL total h2` + pmax(0, pmin(1, h2_obs_polygenic)))

  fwrite(h2_estimates_wide, "h2_estimates_wide.txt", sep="\t", row.names = F, col.names=T)

  h2_estimates_wide <- fread("h2_estimates_wide.txt")

  z_threshold <- qnorm(0.95)
  h2_estimates_wide_ranked <- h2_estimates_wide %>%
    mutate(trans_total_h2 = `Total common variant heritability` - `h2_obs_cis-eQTLs`,
           trans2_total_h2 = `h2_obs_trans-eQTLs` + h2_obs_polygenic,
           trans_total_h2_z = scale(na.omit(trans_total_h2))[,1],
           trans2_total_h2_z = scale(na.omit(trans2_total_h2))[,1],
           cis_total_h2 = `h2_obs_cis-eQTLs`,
           cis_total_h2_z = scale(na.omit(cis_total_h2))[,1],
           cis_larger_than_average = cis_total_h2_z > z_threshold,
           trans2_larger_than_average = trans2_total_h2_z > z_threshold,
           trans_larger_than_average = trans_total_h2_z > z_threshold
    )

  fwrite(h2_estimates_wide_ranked, "h2_estimates_wide_ranked.txt", sep="\t", row.names=F, col.names=T)

  ordered_genes <- h2_estimates %>%
    group_by(phenotype) %>%
    summarise(total_h2 = sum(h2)) %>%
    arrange(total_h2) %>% pull(phenotype)

  h2_estimates_ordered <- h2_estimates %>%
    mutate(gene_id_ranked = as.numeric(ordered(phenotype, levels=ordered_genes)))

  p1 <- ggplot(h2_estimates_ordered, aes(gene_id_ranked, h2, fill = type)) +
    geom_col(position=position_stack()) +
    scale_fill_manual(values = c("polygenic" = "#D55E00", "cis-eQTLs" = "#56B4E9", "trans-eQTLs" = "#F0E442")) +
    xlab("Gene index") + ylab("Common variant heritability") + theme(legend.position = "None")

  p2 <- ggplot(h2_estimates_ordered, aes(h2_obs, fill = type)) +
    geom_histogram() +
    scale_fill_manual(values = c("polygenic" = "#D55E00", "cis-eQTLs" = "#56B4E9", "trans-eQTLs" = "#F0E442")) +
    xlab("Common variant heritability") + ylab("Number of genes") +
    facet_wrap(vars(type), scales="free_x") + theme(legend.position = "top")

  pcomb <- p2 / p1
  pcomb
  ggsave("heritability.pdf", pcomb)

  h2_estimates %>%
    group_by(type) %>%
    summarise(h2_mean = mean(h2_obs), h2_median = median(h2_obs), sum(h2_p < 0.05) / n())

  ldsc %>%
    summarise(h2_mean = mean(h2_obs), h2_median = median(h2_obs), sum(h2_p < 0.05) / n())

  p1 <- ggplot(ldsc, aes(h2_obs, h2_obs_se)) + geom_point()
  p2 <- ggplot(ldsc, aes(h2_obs_lim)) + geom_histogram()

  pggarrange(p1, p2, ncol=1, align="v")


  # Get enrichments of genes with high cis and trans-heritability
  cis_enrichments <- fread("cis-heritability_enrichments.txt") # These files are generated with Toppfun
  trans_enrichments <- fread("trans-heritability_enrichments.txt") # These files are generated with Toppfun

  all_heritability_enrichments <- bind_rows(`cis-heritability` = cis_enrichments, `trans-heritability` = trans_enrichments, .id = "run")

  fwrite(all_heritability_enrichments, "heritability_enrichments.txt", sep = "\t", row.names = F, col.names = T)

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}