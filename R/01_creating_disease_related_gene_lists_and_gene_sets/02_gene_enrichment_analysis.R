library(data.table)
library(tidyverse)
library(enrichR)
library(dplyr)

input_file_path <-
  "data/01_geneLists/01_knownDiseaseRelatedGenes/"

 output_file_path <-
   "data/01_geneLists/02_geneEnrichmentResults/"

combined_gene_list <-
  paste0(input_file_path, "combinedGeneList.csv") %>%
  fread() %>%
  as.data.frame() %>%
  dplyr::select(x) %>%
  as.list() %>%
  unlist() %>%
  unname() %>%
  c()

dbs <- listEnrichrDbs()

enriched <- enrichr(combined_gene_list, dbs$libraryName)

dir.create(
  output_file_path,
  recursive = TRUE,
  showWarnings = FALSE
)

for (n in seq(length(enriched))) {
  write.csv(enriched[[n]],
            paste0(
              output_file_path,
              names(enriched)[n],
              ".csv"
            ))
}
