library(data.table)
library(dplyr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]
mfa <- commandArgs(trailingOnly = TRUE)[2] %>% as.numeric() %>%
  format(scientific = FALSE)

# input_file_name <- "als"
# mfa <- 0.01

input_file_path <-
  paste0("data/07_quality_controlled_gwas_data/",
         input_file_name,
         "/",
         mfa,
         "/quality_controlled_gwas_data.tsv")


output_file_path <-
  paste0("data/12_ldsc_results/",
         input_file_name,
         "/",
         mfa,
         "/")

dir.create(output_file_path, recursive = TRUE, showWarnings = FALSE)

output_file_path <- paste0(output_file_path, "quality_controlled_gwas_data.tsv")

input_file_path %>%
  fread() %>%
  rename(MarkerName = rsid,
         Allele1 = effect_allele,
         Allele2 = other_allele,
         Freq.Allele1.HapMapCEU = effect_allele_frequency,
         p = p_value,
         N = N_effective) %>%
  select(MarkerName, Allele1, Allele2, Freq.Allele1.HapMapCEU, p, N) %>%
  fwrite(output_file_path, sep = "\t")
