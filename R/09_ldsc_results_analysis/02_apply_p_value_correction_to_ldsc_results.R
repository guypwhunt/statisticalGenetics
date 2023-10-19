library(data.table)
library(dplyr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]
mfa <- commandArgs(trailingOnly = TRUE)[2] %>% as.numeric() %>%
  format(scientific = FALSE)

# input_file_name <- "als"
# mfa <- 0.0025 %>%
#   format(scientific = FALSE)

input_file_path <-
  paste0("data/13_consolidated_ldsc_results/",
         input_file_name,
         "/",
         mfa,
         "/")

input_file_names <-
  list.files(input_file_path)

output_file_path <-
  paste0("data/14_p_value_adjusted_ldsc_results/",
         input_file_name,
         "/",
         mfa,
         "/")

dir.create(output_file_path,
           showWarnings = FALSE,
           recursive = TRUE)

for (file in input_file_names) {
  input_file_name <-
    paste0(input_file_path, file)

  gene_set_results <-
    input_file_name %>%
    fread(fill = TRUE) %>%
    arrange(Enrichment_p) %>%
    filter(
      !gene_set %in% c(
        "diseaseAssociatedGenes_diseaseAssociatedGenes",
        "diseaseAssociatedGenes_gwasAndDiseaseAssociatedGenes",
        "diseaseAssociatedGenes_gwasAssociatedGenes"
      )
    ) %>%
    relocate(gene_set)

  gene_set_results$Adjust_P <-
    p.adjust(gene_set_results$Enrichment_p,
             method = "BH")

  fwrite(gene_set_results,
         paste0(output_file_path, file))
}
