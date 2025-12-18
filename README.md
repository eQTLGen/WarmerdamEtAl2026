# Code accompanying the manuscript: "_Trans_-eQTLs reveal the architecture of human gene regulatory networks"

This repository contains code that is used to calculate the numbers and draw the basis for the figures presented in:
Warmerdam et al., _Trans_-eQTLs reveal the architecture of human gene regulatory networks.

## Content
- CompareEqtlsRelevance: Code to run analysis on activity of cis-eQTLs
- ComparisonPhase1: Code to compare our eQTLs to eQTLGen phase 1
- ExtractFinemappedeQTLs: Extract finemapped eQTLs from INTERVAL sQTL data
- MakeCohortOverview: Draw map with locations of cohorts
- MakeFinemappingResults: Analysis for the interpretation of fine-mapped eQTLs
- GwasColocalizations: Script to prepare figures on the GWAS colocalization results
- NetworkAnalysis: Analysis to identify trait-relevant genes linked to multiple GWAS risk alleles through the cis-trans network
- ReplicationAnalyses: Scripts to calculate Rb estimates of our eQTLs in replication datasets

## Running the code

The majority of the code can be run by, defining an environment variable in an `.Renviron` file named `EQTLGEN_RESOURCES`.
This `EQTLGEN_RESOURCES` should be a path to a tab-separated table with variable names and paths to input files.

An example of such a file is `resources/file_list.tsv`