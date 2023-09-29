library(data.table)
library(dplyr)

input_file_path <-
  "data/02_magma/08_geneSetAnalysis/"

input_file_names <-
  list.files(input_file_path, pattern = ".gsa.out")

output_file_path <-
  "data/02_magma/09_pValueAdjustedGeneSetAnalysis/"

for(file in input_file_names) {
  input_file_name <-
    paste0(input_file_path, file)

  gene_set_results <-
    input_file_name %>%
    fread(skip = 4,
          fill = TRUE) %>%
    arrange(P) %>%
    filter(
      !FULL_NAME %in% c(
        "diseaseRelatedGenes_gwasCatalogueGenes.csv",
        "diseaseRelatedGenes_gwasCataloguePlusOrginalGenes.csv",
        "diseaseRelatedGenes_combinedGeneList.csv"
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


