library(data.table)
library(dplyr)

input_file_path <-
  "data/02_magma/08_geneSetAnalysis/"

output_file_path <-
  "data/02_magma/09_pValueAdjustedGeneSetAnalysis/"

fileNames <- c(0, 1, 2.5, 5, 7.5, 10)

gene_set_results <- lapply(fileNames, function(fileName) {
  temp_gene_set_results <-
    paste0(input_file_path, fileName, ".gsa.out") %>%
    fread(skip = 4,
          fill = TRUE) %>%
    arrange(P) %>%
    filter(
      !FULL_NAME %in% c(
        "diseaseRelatedGenes_gwasCatalogueGenes.csv",
        "diseaseRelatedGenes_gwasCataloguePlusOrginalGenes.csv",
        "diseaseRelatedGenes_combinedGeneList.csv"
      )
    )

  temp_gene_set_results$Adjust_P <-
    p.adjust(temp_gene_set_results$P,
             method = "BH")

  return(temp_gene_set_results)
})

names(gene_set_results) <- fileNames

dir.create(output_file_path, showWarnings = FALSE, recursive = TRUE)
