library(data.table)
library(biomaRt)
library(dplyr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

# input_file_name <- "als"

input_file_path <-
  paste0(
    "data/02_gene_set_enrichment/04_gene_enrichment_term_genes/",
    input_file_name,
    "/"
  )

output_file_path <-
  paste0("data/06_format_gene_sets/",
         input_file_name,
         "/01_gene_sets/")


dglinker_genes_input_file_path <-
  paste0("data/03_dglinker_results_analysis/02_dglinker_gene_list/",
         input_file_name,
         "/")

additional_genes_input_file_path <-
  paste0("data/04_additional_gene_sets/", input_file_name, "/")

directories <- sub('.*\\/', '', list.dirs(input_file_path))

directories <- directories[directories != ""]

ensembl <-
  useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

for (directory in directories) {
  input_directory <- paste0(input_file_path, directory, "/")

  output_directory <- paste0(output_file_path, directory, "/")

  dir.create(output_directory,
             showWarnings = FALSE,
             recursive = TRUE)

  input_files <- list.files(input_directory)

  for (input_file in input_files) {
    temp_input_file_path <- paste0(input_directory, input_file)

    gene_list <- fread(temp_input_file_path, header = FALSE) %>%
      as.list() %>%
      unname() %>%
      unlist()

    gene_list <- getBM(
      attributes = c('entrezgene_id'),
      filters = 'hgnc_symbol',
      values = gene_list,
      mart = ensembl
    ) %>%
      unique()

    fwrite(
      gene_list,
      paste0(output_directory,
             input_file),
      row.names = FALSE,
      col.names = FALSE
    )

  }
}

if (file.exists(dglinker_genes_input_file_path)) {
  dglinker_file_names <- list.files(dglinker_genes_input_file_path)

  for (dglinker_file_name in dglinker_file_names) {
    dglinker_file_path <- paste0(dglinker_genes_input_file_path,
                                 dglinker_file_name)

    dglinker_data <-
      fread(dglinker_file_path, header = FALSE) %>% as.list() %>% unname() %>% unlist()

    output_directory <- paste0(output_file_path, "dglinker/")

    dir.create(output_directory,
               showWarnings = FALSE,
               recursive = TRUE)

    dglinker_data <- getBM(
      attributes = c('entrezgene_id'),
      filters = 'hgnc_symbol',
      values = dglinker_data,
      mart = ensembl
    ) %>%
      unique()

    fwrite(
      dglinker_data,
      paste0(output_directory,
             dglinker_file_name),
      row.names = FALSE,
      col.names = FALSE
    )
  }
}

if (file.exists(additional_genes_input_file_path)) {
  additional_genes_file_names <-
    list.files(additional_genes_input_file_path)

  for (additional_genes_file_name in additional_genes_file_names) {
    additional_genes_file_path <-
      paste0(additional_genes_input_file_path,
             additional_genes_file_name)

    additional_genes_data <-
      fread(additional_genes_file_path, header = FALSE) %>% as.list() %>% unname() %>% unlist()

    output_directory <- paste0(output_file_path, "additionalGenes/")

    dir.create(output_directory,
               showWarnings = FALSE,
               recursive = TRUE)

    additional_genes_data <- getBM(
      attributes = c('entrezgene_id'),
      filters = 'hgnc_symbol',
      values = additional_genes_data,
      mart = ensembl
    ) %>%
      unique()

    fwrite(
      additional_genes_data,
      paste0(output_directory,
             additional_genes_file_name),
      row.names = FALSE,
      col.names = FALSE
    )
  }
}
