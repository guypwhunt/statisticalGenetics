library(data.table)
library(tidyverse)
library(enrichR)
library(dplyr)
library(tools)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

# input_file_name <- "als.csv"

input_file_path <-
  "data/01_data_input/01_disease_related_genes/"

output_file_path <- paste0(
    "data/02_gene_set_enrichment/01_gene_enrichment_results/",
    tools::file_path_sans_ext(input_file_name) , "/"
)

combined_gene_list <-
  paste0(input_file_path, input_file_name) %>%
  fread() %>%
  as.data.frame() %>%
  dplyr::select(x) %>%
  as.list() %>%
  unlist() %>%
  unname() %>%
  c()

dbs <- listEnrichrDbs()

enriched <-
  enrichr(
    combined_gene_list,
    c(
      "GO_Biological_Process_2023",
      "GO_Cellular_Component_2023",
      "GO_Molecular_Function_2023"
    )
  )

dir.create(output_file_path,
           recursive = TRUE,
           showWarnings = FALSE)

for (n in seq(length(enriched))) {
  write.csv(
    enriched[[n]],
    paste0(
      output_file_path,
      names(enriched)[n],
      ".csv"
    )
  )
}
