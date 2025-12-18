#!/usr/bin/env Rscript

# Load libraries
library(data.table)
library(rtracklayer)
library(argparse)
library(stringr)
library(tidyverse)
library(dtplyr)
library(arrow)

# Declare constants
SOFTWARE_VERSION <- 2.0

EQTLGEN_P1_FORMAT_INFO <- list(
  "gene" = "Gene",
  "snp" = "SNP",
  "chr" = "SNPChr",
  "bp" = "SNPPos",
  "a1" = "AssessedAllele",
  "a2" = "OtherAllele",
  "effect" = "Zscore",
  "p" = "Pvalue",
  "n" = "NrSamples",
  "eaf" = "AlleleB_all",
  "build" = "hg19")

# Declare function definitions

# Define a function to check if required columns exist and are non-NA
check_columns <- function(data, required_columns) {
  missing_columns <- setdiff(required_columns, names(data))  # Find missing columns
  names(missing_columns) <- names(required_columns)[required_columns %in% missing_columns] # Add names to missing columns

  # Check for missing columns
  if (length(missing_columns) > 0) {
    message(sprintf("The following required columns are missing: ,", paste(names(missing_columns), collapse = ", ")))
  } else {
    message("All required columns are present.")
  }

  columns_with_na <- required_columns[sapply(data[, ..required_columns], function(col) any(is.na(col)))]  # Columns containing NA
  columns_with_all_na <- required_columns[sapply(data[, ..required_columns], function(col) all(is.na(col)))]  # Columns containing NA

  # Check for columns with NA values
  if (length(columns_with_na) > 0) {
    message("The following columns contain some NA values: ", paste(columns_with_na, collapse = ", "))
  } else {
    message("All required columns have non-NA values.")
  }

  # Check for columns with NA values
  if (length(columns_with_all_na) > 0) {
    message("The following columns contain just NA values: ", paste(columns_with_all_na, collapse = ", "))
  }

  # If everything is good
  if (length(missing_columns) == 0 && length(columns_with_all_na) == 0) {
    message("All checks passed. Data is ready for further processing.")
    return(TRUE)
  } else {
    warning("Some issues detected with the data. Please review the messages above.")
    print(head(data))
    return(FALSE)
  }
}

# See Zhang et al, Nature Genetics volume 48, pages481–487(2016) for details
# Formula taken from Supplementary Information
z_to_beta <- function(z, maf, n) {
  chi <- z * z
  a <- 2.0 * maf * (1.0 - maf) * (n + chi)
  beta <- z / sqrt(a)
  se <- 1.0 / sqrt(a)

  return(list("beta" = beta, "se" = se))
}

# Harmonize each GWAS file
harmonize_eqtlgen_p1_files <- function(gwas_data, format_info, liftover_chain = NULL, output_z_scores = FALSE, output_sample_size = TRUE, sample_size = NULL) {
  build <- format_info$build

  column_types <- list(numeric = c(format_info$p),
                       integer = c(format_info$bp),
                       character = c(format_info$chr, format_info$snp, format_info$a1, format_info$a2))

  # Load the summary statistics file using data.table
  original_length <- nrow(gwas_data)

  gwas_data[, (column_types[["numeric"]]) := lapply(.SD, as.numeric), .SDcols = column_types[["numeric"]]]
  gwas_data[, (column_types[["integer"]]) := lapply(.SD, as.numeric), .SDcols = column_types[["integer"]]]

  required_columns <- c("chr" = format_info$chr, "bp" = format_info$bp, "snp" = format_info$snp,
                        "a1" = format_info$a1, "a2" = format_info$a2, "p" = format_info$p, "n" = format_info$n,
                        "eaf" = format_info$eaf)

  if (format_info$snp == "" | !(format_info$snp %in% names(gwas_data))) {
    gwas_data[, "snp" := paste(get(required_columns['chr']),
                               get(required_columns['bp']),
                               get(required_columns['a1']),
                               get(required_columns['a2']), sep="_")]
    required_columns["snp"] <- "snp"
  }

  gwas_data[, required_columns['chr'] := gsub("^0+", "", get(required_columns['chr']))]

  check_columns(gwas_data, required_columns)

  # Rename columns to a standardized format
  setnames(gwas_data, old = required_columns,
           new = c("CHR", "BP", "SNP", "A1", "A2", "p", "N", "eaf"))

  # Check and convert p-values
  if (max(gwas_data$p, na.rm=T) > 1 & str_detect(tolower(format_info$p), "log")) {
    message("Detected -log10 p-values. Converting...")
    gwas_data[, p := 10^(-gwas_data$p)]
  } else if (max(gwas_data$p, na.rm=T) > 1) {
    warning("Max p-value is greater than 1, but column name does not indicate log transformation. Check this!")
  } else if (str_detect(tolower(format_info$p), "log")) {
    warning("Max p-value is not greater than 1, but column name indicates log transformation. Check this!")
  }

  if (min(gwas_data$p, na.rm=T) <= 0) {
    message(sprintf("Detected %s variants with P <= 0", sum(gwas_data$p <= 0, na.rm=T)))
    gwas_data[p <= 0, p := .Machine$double.xmin * .Machine$double.eps * 2]
  }

  beta_se <- z_to_beta(gwas_data$Zscore, gwas_data$eaf, gwas_data$N)
  gwas_data[, EFFECT := beta_se$beta]
  gwas_data[, se := beta_se$se]

  columns_to_check <- c("CHR", "BP", "SNP", "A1", "A2", "p", "se", "EFFECT", "N", "Zscore")

  print(gwas_data)

  gwas_data <- gwas_data[rowSums(is.na(gwas_data[, columns_to_check, with = FALSE])) == 0, ]
  removed_na_length <- nrow(gwas_data)

  # Harmonize to build 38 if necessary
  if (build != "hg38" && !is.null(liftover_chain)) {
    gwas_data <- lift_to_build_38(gwas_data, liftover_chain)
  }

  lifted_length <- nrow(gwas_data)

  return(list(gwas_dt = gwas_data,
              original_length = original_length,
              removed_na_length = removed_na_length,
              lifted_length = lifted_length))
}


# Function to lift genome build to hg38 using rtracklayer (LiftOver)
lift_to_build_38 <- function(gwas_data, liftover_chain_file) {
  # Create a GRanges object for liftover
  gwas_data[, index := seq_along(gwas_data$SNP)]

  gr <- GRanges(
    seqnames = Rle(paste0("chr", str_remove(gwas_data$CHR, "^chr"))),
    ranges = IRanges(start = gwas_data$BP, end = gwas_data$BP),
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
  gwas_data_filtered[, CHR := sub("chr", "", new_coords_dt[["seqnames"]])]
  gwas_data_filtered[, BP := new_coords_dt[["start"]]]
  gwas_data_filtered[ , index := NULL, ]

  return(gwas_data_filtered)
}


check_if_skip <- function(log_file, format_info, gwas_file) {
  if (file.exists(log_file)) {
    log <- fread(log_file)
    print(log)
    comparison_results <- c(log$software_version == SOFTWARE_VERSION, log$original_file == gwas_file,
                            log$format_info_skip_lines == format_info$skip_lines,
                            log$format_info_chromosome == format_info$chromosome,
                            log$format_info_bp == format_info$bp,
                            log$format_info_a1 == format_info$a1,
                            log$format_info_a2 == format_info$a2,
                            log$format_info_snp == format_info$snp,
                            log$format_info_p == format_info$p,
                            log$format_info_z_score == format_info$z_score,
                            log$format_info_effect == format_info$effect,
                            log$format_info_se == format_info$se)
    if (
      all(c("software_version", "format_info_chromosome", "format_info_bp", "format_info_a1",
            "format_info_a2","format_info_snp","format_info_p","format_info_effect","format_info_se") %in% names(log))
        & all(comparison_results, na.rm=T)) {
      return(1)
    }
    return(-1)
  }
  return(0)
}


# Main function to run the pipeline
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = TRUE)
  }

  # Set up argument parser
  parser <- ArgumentParser(description = "Harmonize GWAS summary statistics files.")

  # Add arguments
  parser$add_argument("--out-dir", type = "character", required = TRUE,
                      help = "Base directory where the processed summary statistics files will be stored.")
  parser$add_argument('--gwas-file', type = 'character', nargs="+",
                      help = 'A tab-delimited file with GWAS summary statistics.')
  parser$add_argument('--trait-id', type = 'character', help = "trait id, to include in output file")
  parser$add_argument("--af-file", type = "character", required = TRUE,
                      help = 'Allele frequency file')
  parser$add_argument("--liftover-chain", type = "character", required = TRUE,
                      help = "Path to the chain file for liftover (e.g., hg19ToHg38.over.chain.gz).")
  parser$add_argument("--sample-size", "-n", type = "integer", required = FALSE,
                      help = "Sample size of the gwas (is this is not in the GWAS data itself)")
  parser$add_argument("--eQTLGen-index", type = "character", required = TRUE,
                      help = "extended eQTLGen parquet index file (with unique mapping ids).")

  # Parse arguments
  args <- parser$parse_args(argv)

  eqtlgen_index_file <- args$eQTLGen_index
  eqtl_paths <- args$gwas_file
  af_path <- args$af_file
  liftover_chain_file <- args$liftover_chain
  out_dir <- args$out_dir
  trait_id <- args$trait_id

  # Get format information for the current GWAS file
  format_info <- EQTLGEN_P1_FORMAT_INFO

  gwas_data_original <- bind_rows(mapply(function(gwas_path) {
    message(sprintf("Reading file: '%s'", gwas_path))
    fread(gwas_path, colClasses = list(character = "GeneChr"))
  }, eqtl_paths, SIMPLIFY=F, USE.NAMES = F))

  eqtlgen_index_extended <- arrow::read_parquet(eqtlgen_index_file)
  setkey(eqtlgen_index_extended, mapping_id)

  allele_frequencies <- fread(af_path, select = c("SNP", "AlleleB_all"))

  gwas_data <- copy(gwas_data_original)

  gwas_data <- allele_frequencies[gwas_data,
           on = "SNP"]

  # Harmonize the file
  harmonized_data <- harmonize_eqtlgen_p1_files(gwas_data, format_info, liftover_chain_file, output_z_scores = FALSE, output_sample_size = TRUE, sample_size = sample_size)
  harmonized_gwas_data <- harmonized_data$gwas_dt
  original_length <- harmonized_data$original_length
  removed_na_length <- harmonized_data$removed_na_length
  lifted_length <- harmonized_data$lifted_length

  message(sprintf("Original rows: %i\nAfter removing rows where some of the remapped columns have NA's: %i\nAfter liftOver: %i\n", original_length, removed_na_length, lifted_length))

  # Save the harmonized data to a new file
  output_file <- file.path(out_dir, sprintf("%s.sumstats.harmonized.gz", trait_id))
  fwrite(harmonized_gwas_data[, .(Gene, SNP, CHR, BP, A1, A2, Zscore, EFFECT, p, se, N)], output_file, sep = "\t", row.names=F, col.names=T, quote=F)
  message(sprintf("Harmonized %s -> %s", eqtl_paths, output_file))
  message("Starting merging...")

  harmonized_gwas_data[, mapping_id := paste0(CHR, "_", BP, "_", pmin(A1, A2), "_", pmax(A1, A2))]
  harmonized_data_subset <- harmonized_gwas_data[, .(Gene, mapping_id, A1, A2, Zscore, EFFECT, p, se, eaf, N)]

  setkey(harmonized_data_subset, mapping_id)
  merged_gwas_data <- merge(harmonized_data_subset,eqtlgen_index_extended, all=F)
  merged_length <- nrow(merged_gwas_data)
  message(sprintf("After merging: %i", merged_length))

  merged_gwas_data[, beta := merged_gwas_data[["EFFECT"]]]
  merged_gwas_data[, z := merged_gwas_data[["Zscore"]]]
  merged_gwas_data[eff_allele != A1,"beta"] <- merged_gwas_data[eff_allele != A1,beta*-1]
  merged_gwas_data[eff_allele != A1,"eaf"] <- merged_gwas_data[eff_allele != A1,eaf*-1]
  merged_gwas_data[eff_allele != A1,"z"] <- merged_gwas_data[eff_allele != A1,Zscore*-1]

  swapped_rows <- nrow(merged_gwas_data[beta != EFFECT])
  message(sprintf("Swapped beta sign: %i / %i (%.3f%%)", swapped_rows, merged_length, swapped_rows/merged_length*100))

  merged_gwas_data[, c("mapping_id", "EFFECT", "Zscore") := NULL, ]

  merged_gwas_data[, gwas_id := trait_id]

  write_dataset(merged_gwas_data, file.path(out_dir, "eqtlgen_parquet_dataset"),
                format="parquet", partitioning="gwas_id", existing_data_behavior="delete_matching")

  log <- tibble(trait_id = args$pqtl_id,
                original_file = args$pqtl_folder,
                format = "eQTLGen_phase1",
                original_length = original_length,
                removed_na_length = removed_na_length,
                lifted_length = lifted_length,
                merged_length = merged_length,
                swapped_rows = swapped_rows,
                format_info = list(format_info),
                software_version = SOFTWARE_VERSION
  ) %>% unnest_wider(format_info, names_sep="_")

  log_file <- file.path(args$out_dir, sprintf("%s.log.txt", trait_id))

  fwrite(log, log_file)
}

# Entry point check
if (sys.nframe() == 0 && !interactive()) {
  main()
}
