#!/bin/bash
#SBATCH --job-name=gene_prioritization
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=96gb
#SBATCH --output=gene_prioritization_%j.out
#SBATCH --error=gene_prioritization_%j.err

# ---- Absolute file paths based on the original R script ----

ml R

CONFIG="resources.txt"

OUTPUT_DIR="network_analyses"

HGNC="hgnc_complete_set.txt"
CYTOBAND="cytoBand.txt.gz"
CIS_TRANS="" # Cis to trans causal relationships
GWAS_DATASET="" # GWAS dataset in parquet format
UNCORRELATED_VARIANTS="" # A set of independent variants
GWAS_EQTL_COLOC="" # eQTL to GWAS colocalization output
EQTL_MATRIX="" # A matrix with eQTL Z-scores for all unique fine-mapped variants and all genes
SC_REPLICATION="" # Folder with single-cell eQTLGen replication summary statistics

OUT_PREFIX="" # Output prefix

PLOT_GENE="ABCA1"
PLOT_GWAS="Low_Density_Lipoprotein_Cholesterol_levels"

cd $OUTPUT_DIR

# ---- Run the pipeline with config + all overrides ----

Rscript network_analyses_coloc.R \
  --config "$CONFIG" \
  --hgnc_file "$HGNC" \
  --cytoband_file "$CYTOBAND" \
  --cis_trans_file "$CIS_TRANS" \
  --gwas_dataset "$GWAS_DATASET" \
  --uncorrelated_variant_path "$UNCORRELATED_VARIANTS" \
  --gwas_eqtl_coloc "$GWAS_EQTL_COLOC" \
  --eqtl_matrix "$EQTL_MATRIX" \
  --sc_replication "$SC_REPLICATION" \
  --disable_sc_cond \
  --output_prefix $OUT_PREFIX \
  --focal_gene "$PLOT_GENE" \
  --focal_gwas "$PLOT_GWAS"
