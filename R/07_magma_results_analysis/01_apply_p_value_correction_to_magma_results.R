library(data.table)
library(dplyr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]
mfa <- commandArgs(trailingOnly = TRUE)[2] %>% as.numeric() %>%
  format(scientific = FALSE)

# input_file_name <- "als"

input_file_path <-
  paste0("data/09_magma_results/", input_file_name, "/",
         mfa,
         "/")

input_file_names <-
  list.files(input_file_path, pattern = ".gsa.out")

output_file_path <-
  paste0("data/10_p_value_adjusted_magma_results/",
         input_file_name,
         "/",
         mfa,
         "/")

for (file in input_file_names) {
  input_file_name <-
    paste0(input_file_path, file)

  gene_set_results <-
    input_file_name %>%
    fread(skip = 4,
          fill = TRUE) %>%
    arrange(P) %>%
    filter(
      !FULL_NAME %in% c(
        "diseaseAssociatedGenes_diseaseAssociatedGenes.csv",
        "diseaseAssociatedGenes_gwasAndDiseaseAssociatedGenes.csv",
        "diseaseAssociatedGenes_gwasAssociatedGenes.csv"
      )
    ) %>%
    relocate(FULL_NAME) %>%
    dplyr::select(!VARIABLE)

  gene_set_results$Adjust_P <-
    p.adjust(gene_set_results$P,
             method = "BH")

  dir.create(output_file_path,
             showWarnings = FALSE,
             recursive = TRUE)

  fwrite(gene_set_results,
         paste0(output_file_path, file, ".csv"))
}
