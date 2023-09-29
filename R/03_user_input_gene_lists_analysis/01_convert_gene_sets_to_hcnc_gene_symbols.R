library(data.table)
library(tidyverse)
library(dplyr)
library(stringr)
library(readxl)
library(biomaRt)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

# input_file_name <- "als"

input_file_path <-
  paste0("data/01_data_input/04_additional_gene_sets/",
         input_file_name,
         "/")

output_file_path <-
  paste0("data/04_additional_gene_sets/",
         input_file_name,
         "/")

dir.create(output_file_path, showWarnings = FALSE, recursive = TRUE)

file_names <- list.files(input_file_path)

if (length(file_names) > 0) {
  ensembl <-
    useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

  for (file_name in file_names) {
    gene_set <-
      fread(paste0(input_file_path, file_name), header = FALSE, fill=TRUE) %>%
      as.list() %>%
      unlist()

    gene_set <- getBM(
      attributes = c('hgnc_symbol'),
      filters = 'hgnc_symbol',
      values = gene_set,
      mart = ensembl
    ) %>%
      as.list() %>%
      unlist() %>%
      unique() %>%
      as.data.frame()

    fwrite(gene_set, paste0(output_file_path, file_name), row.names = FALSE, col.names = FALSE)
  }
}
