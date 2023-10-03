library(biomaRt)
library(dplyr)
library(tidyverse)
library(data.table)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

# input_file_name <- "als"

input_file_path <-
  paste0(
    "data/02_gene_set_enrichment/03_significant_non_redundant_gene_enrichmentResults/",
    input_file_name,
    "/"
  )

output_file_path <-
  paste0(
    "data/02_gene_set_enrichment/04_gene_enrichment_term_genes/",
    input_file_name,
    "/"
  )

dir.create(output_file_path,
           showWarnings = FALSE,
           recursive = TRUE)


go_file_names <-
  list.files(
    path = input_file_path,
    pattern = "GO_",
    all.files = FALSE,
    full.names = FALSE
  ) %>% rev ()

for (go_file_name in go_file_names) {
  go_file_paths <- paste0(input_file_path, go_file_name)

  go_data_frame <- read.csv(go_file_paths, header = TRUE)

  go_data_frame$GO_DB <-
    str_replace(go_file_name, ".csv", "")

  ensembl <-
    useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  go_ids <- unique(go_data_frame$ID)

  go_terms <-
    str_replace_all(unique(go_data_frame$Term), "[[:punct:]]", "")

  go_db <- unique(go_data_frame$GO_DB)

  gene_data_frame <-
    lapply(go_ids, function(i) {
      getBM(
        attributes = 'hgnc_symbol',
        filters = "go_parent_term",
        uniqueRows = TRUE,
        values = i,
        mart = ensembl
      ) %>%
        as.list() %>%
        unlist() %>%
        unique() %>%
        as.data.frame()
    })

  dir.create(paste0(output_file_path, go_db),
             showWarnings = FALSE,
             recursive = TRUE)

  for (gene_data_frame_number in seq(length(gene_data_frame))) {
    temp_gene_data_frame <- gene_data_frame[[gene_data_frame_number]]

    fwrite(
      temp_gene_data_frame,
      paste0(output_file_path, go_db, "/", go_terms[gene_data_frame_number], ".csv"),
      row.names = FALSE, col.names = FALSE
    )
  }
  message(go_db)
}
