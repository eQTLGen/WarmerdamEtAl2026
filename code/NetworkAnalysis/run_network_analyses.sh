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

CONFIG="/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/users/umcg-rwarmerdam/resources.txt"

OUTPUT_DIR="/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/freeze3/Interpretation/network_analyses"

HGNC="/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/freeze3/Interpretation/network_analyses/hgnc_complete_set.txt"
CYTOBAND="/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/freeze3/Interpretation/network_analyses/cytoBand.txt.gz"
CIS_TRANS="/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/freeze3/Interpretation/network_analyses/CausalCisToTransGenesBasedOnFinemappingResultsCisAndTransSNPNamesBasedOnColocalization_20250722.txt"
INDEP_VARIANTS="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/users/umcg-rwarmerdam/independent_variants.txt"
GWAS_DATASET="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/analyses/gwas_processing/gwas_data/eqtlgen_parquet_dataset"
UNCORRELATED_VARIANTS="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/extractions/uncorrelated_variants_2025_01_22/merged.prune.in"
GWAS_EQTL_COLOC="eqtl_gwas_coloc_cleaned_2025-07-30.txt"
EQTL_MATRIX="/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/freeze3/Interpretation/extractions/significant_matrix_finemapped_all_20250721/exported_matrix/mat.z.txt.gz"
SC_REPLICATION='/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/freeze3/Interpretation/eqltVariantAnalysis/singleCellReplications/Meta_v5/UMCG_v2-UMCG_v3-oneK1k-Multiome-Wijst.*.incPerm.tsv'

GRAPH_FILE="directed_gene_network_with_disease_associations_as_igraph_coloc_AllEqtls.rds"
DIST_MATRIX="dist_matrix_coloc_final_cis_trans_AllEqtls.rds"
PATH_TABLE="path_table_coloc_final_cis_trans_AllEqtls.rds"
PATH_TABLE_FILTERED="paths_filtered_coloc_AllEqtls.rds"
OUTPUT_PERMUTATION="gene_prioritization_1000_random_variants_permutations_coloc_AllEqtls_20250901.txt.gz"

PLOT_GENE="ABCA1"
PLOT_GWAS="Low_Density_Lipoprotein_Cholesterol_levels"

cd $OUTPUT_DIR
export R_LIBS="/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/users/umcg-rwarmerdam/R/x86_64-pc-linux-gnu-library/4.4"

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
  --graph_file "$GRAPH_FILE" \
  --dist_matrix_file "$DIST_MATRIX" \
  --path_table_file "$PATH_TABLE" \
  --output_path_table_filtered "$PATH_TABLE_FILTERED" \
  --output_permutation_results "$OUTPUT_PERMUTATION" \
  --focal_gene "$PLOT_GENE" \
  --focal_gwas "$PLOT_GWAS"
