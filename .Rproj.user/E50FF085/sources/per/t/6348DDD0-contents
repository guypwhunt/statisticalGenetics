library(biomaRt)
library(dplyr)
library(tidyverse)

dir.create("data/geneLists/finalGeneLists",
           showWarnings = FALSE,
           recursive = TRUE)

file_path <- "data/geneLists/filteredEnrichmentResults"

go_file_names <-
  list.files(
    path = file_path,
    pattern = "GO_",
    all.files = FALSE,
    full.names = FALSE
  ) %>% rev ()

for (go_file_name in go_file_names) {
  go_file_paths <- paste0(file_path, "/", go_file_name)

  go_data_frame <- read.csv(go_file_paths, header = TRUE)

  go_data_frame$GO_DB <-
    str_replace(go_file_name, ".csv", "")

  ensembl <-
    useMart("ensembl", dataset = "hsapiens_gene_ensembl") #uses human ensembl annotations

  go_ids <- unique(go_data_frame$ID)
  go_terms <-
    str_replace_all(unique(go_data_frame$Term), "[[:punct:]]", "")
  go_db <- unique(go_data_frame$GO_DB)

  gene_data_frame <-
    lapply(go_ids, function(i) {
      getBM(
        attributes = c('hgnc_symbol'#, 'chromosome_name', 'start_position', 'end_position'
                       ),
        # filters = 'go',
        filters = "go_parent_term",
        uniqueRows = TRUE,
        values = i,
        mart = ensembl
      )
    })

  dir.create(paste0("data/geneLists/finalGeneLists/", go_db),
             showWarnings = FALSE,
             recursive = TRUE)

  for (gene_data_frame_number in seq(length(gene_data_frame))) {
    temp_gene_data_frame <- gene_data_frame[[gene_data_frame_number]]

    write.csv(
      temp_gene_data_frame,
      paste0("data/geneLists/finalGeneLists/", go_db, "/", go_terms[gene_data_frame_number], ".csv"),
      row.names = FALSE
    )
  }
  message(go_db)
}
