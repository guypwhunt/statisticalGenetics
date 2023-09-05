library(data.table)
library(dplyr)

input_file_path <-
  "data/02_magma/08_geneSetAnalysis"

output_file_path <-
  "data/02_magma/09_pValueAdjustedGeneSetAnalysis/"

fileNames <- c(0, 1, 2.5, 5, 7.5, 10)
fileName <- fileNames[1]

directories <- list.dirs(input_file_path)[2:4]
directory <- directories[1]

for (directory in directories) {
  for (fileName in fileNames) {
    gene_set_results <-
      paste0(directory, "/", fileName, ".gsa.out") %>%
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

    full_output_file_path <-
      paste0(output_file_path, sub('.*\\/', '', directory), "/")


    dir.create(full_output_file_path,
               showWarnings = FALSE,
               recursive = TRUE)

    fwrite(gene_set_results,
           paste0(full_output_file_path, fileName, ".csv"))
  }
}
