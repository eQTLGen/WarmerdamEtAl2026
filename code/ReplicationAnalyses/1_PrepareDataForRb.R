library(data.table)
library(tidyverse)
library(IGUtilityPackage)

variables <- deframe(read_delim(Sys.getenv("EQTLGEN_RESOURCES"), delim = "\t"))

rep_table <- fread(variables["replication_data"])

rep_table$metaLclTstat <- rep_table$metaLclBeta / rep_table$metaLclSe
rep_table$metaLclFdr <- p.adjust(10^-(rep_table$metaLclNlog10p), method = "fdr")

rep_table <- rep_table[, c(1, 21, 46, 3, 4, 70:124, 131:338, 576, 577, 573, 574), with = FALSE]
colnames(rep_table)[1] <- "eQTL"
colnames(rep_table)[c(269, 270, 271, 272)] <- c("meta_LCL_tstat", "meta_LCL_fdr", "meta_LCL_beta", "meta_LCL_se")

beta_columns <- colnames(rep_table)[str_detect(colnames(rep_table), pattern = "beta") | str_detect(colnames(rep_table), pattern = "Beta")][-c(1, 3)]

rep_betas <- rep_table[, colnames(rep_table) %in% c("eQTL", beta_columns), with = FALSE]
rep_betas <- melt(rep_betas, id.vars = "eQTL")

se_columns <- colnames(rep_table)[str_detect(colnames(rep_table), pattern = "_se") | str_detect(colnames(rep_table), pattern = "metaBrainBetaSe")]

rep_ses <- rep_table[, colnames(rep_table) %in% c("eQTL", se_columns), with = FALSE]
rep_ses <- melt(rep_ses, id.vars = "eQTL")

rep_discovery <- rep_table[, c(1, 2, 3, 4, 5), with = FALSE]
colnames(rep_discovery)[5] <- "se"

rep_betas$variable <- str_remove(rep_betas$variable, "Beta")
rep_betas$variable <- str_remove(rep_betas$variable, "_beta")
rep_betas$variable <- str_replace(rep_betas$variable, "metaBrainBeta", "metaBrain")

rep_ses$variable <- str_remove(rep_ses$variable, "BetaSe")
rep_ses$variable <- str_remove(rep_ses$variable, "_se")
rep_ses$variable <- str_replace(rep_ses$variable, "metaBrainSe", "metaBrain")

table(rep_ses$variable %in% rep_betas$variable)
table(rep_betas$variable %in% rep_ses$variable)

rep <- merge(rep_discovery, rep_betas, by = "eQTL", all.x = TRUE)
colnames(rep)[ncol(rep)] <- "rep_beta"
rep <- merge(rep, rep_ses, by = c("eQTL", "variable"), all.x = TRUE)
colnames(rep)[ncol(rep)] <- "rep_se"

# Remove low-quality PER dataset
rep <- rep[variable != "scBrain_PER"]

# Add average effect over sc-RNAseq blood datasets
rep_sc_blood <- rep[str_detect(variable, "scMeta")]
rep_sc_blood$repZ <- rep_sc_blood$rep_beta / rep_sc_blood$rep_se

rep_sc_blood <- rep_sc_blood %>%
  group_by(eQTL) %>%
  summarise(meanZ = mean(repZ, na.rm = TRUE), nr_not_na = length(repZ[!is.na(repZ)]), eaf = unique(eaf)) %>%
  as.data.table()

rep_sc_blood <- rep_sc_blood %>% mutate(eQTL2 = eQTL) %>% separate(eQTL2, into = c("variant", "gene"), sep = "-") %>% as.data.table()

rep_sc_blood$N <- 1323

ConvertedBetas <- deriveBetaSE(rep_sc_blood, SnpCol = "variant", GeneCol = "gene", ZCol = "meanZ", NCol = "N", MafCol = "eaf")
ConvertedBetas$eQTL <- paste(ConvertedBetas$SNP, ConvertedBetas$Gene, sep = "-")
ConvertedBetas$variable <- "scMeta_average"
colnames(ConvertedBetas)[c(3, 4)] <- c("rep_beta", "rep_se")


rep_f <- unique(rep[, c(1, 3, 4, 5, 6), with = FALSE])

ConvertedBetas <- merge(rep_f, ConvertedBetas[, c(5, 6, 3, 4), with = FALSE])
ConvertedBetas <- ConvertedBetas[, c(1, 6, 2:5, 7, 8), with = FALSE]

rep <- rbind(rep, ConvertedBetas)


# Add average effect over sc-RNAseq brain datasets
rep_sc_brain <- rep[str_detect(variable, "scBrain")]

rep_sc_brain$repZ <- rep_sc_brain$rep_beta / rep_sc_brain$rep_se

rep_sc_brain <- rep_sc_brain %>%
  group_by(eQTL) %>%
  summarise(meanZ = mean(repZ, na.rm = TRUE), nr_not_na = length(repZ[!is.na(repZ)]), eaf = unique(eaf)) %>%
  as.data.table()

rep_sc_brain <- rep_sc_brain %>% mutate(eQTL2 = eQTL) %>% separate(eQTL2, into = c("variant", "gene"), sep = "-") %>% as.data.table()

rep_sc_brain$N <- 757

ConvertedBetas <- deriveBetaSE(rep_sc_brain, SnpCol = "variant", GeneCol = "gene", ZCol = "meanZ", NCol = "N", MafCol = "eaf")
ConvertedBetas$eQTL <- paste(ConvertedBetas$SNP, ConvertedBetas$Gene, sep = "-")
ConvertedBetas$variable <- "scBrain_average"
colnames(ConvertedBetas)[c(3, 4)] <- c("rep_beta", "rep_se")

rep_f <- unique(rep[, c(1, 3, 4, 5, 6), with = FALSE])

ConvertedBetas <- merge(rep_f, ConvertedBetas[, c(5, 6, 3, 4), with = FALSE])
ConvertedBetas <- ConvertedBetas[, c(1, 6, 2:5, 7, 8), with = FALSE]

rep <- rbind(rep, ConvertedBetas)

rep <- rep[, c(1, 3, 2, 4:ncol(rep)), with = FALSE]
fwrite(rep, "rb_results/replications_formatted_for_rb_2025-12-17.txt")

# Keep only effects that are independent from each other (gene-based LD-pruning/clumping)

clumped_eQTLs <- fread(variables["clumped_finemapping_results"])
clumped_eQTLs$comb <- paste(clumped_eQTLs$variant, clumped_eQTLs$phenotype, sep = "-")

rep_f <- rep
rep_f$gene <- str_remove(rep_f$eQTL, ".*-")

rep_f <- rep_f[eQTL %in% clumped_eQTLs$comb]

fwrite(rep_f, "rb_results/replications_formatted_for_rb_ClumpedVariants_2025-12-17.txt")
