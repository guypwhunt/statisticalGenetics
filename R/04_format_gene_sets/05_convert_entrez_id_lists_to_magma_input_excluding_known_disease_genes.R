library(data.table)
library(dplyr)
library(stringr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

input_file_name <- "als"

input_file_path <-
  paste0("data/06_format_gene_sets/",
         input_file_name,
         "/01_gene_sets/")

output_file_path <-
  paste0("data/06_format_gene_sets/",
         input_file_name,
         "/02_magmaInput/")

dir.create(output_file_path, showWarnings = FALSE)

directories <- sub('.*\\/', '', list.dirs(input_file_path))

directories <- directories[directories != ""]

magma_input <- list()

genesToRemove <-
  paste0(
    "data/06_format_gene_sets/",
    input_file_name,
    "/01_gene_sets/diseaseAssociatedGenes/gwasAndDiseaseAssociatedGenes.csv"
  ) %>%
  fread(header = FALSE) %>%
  as.list() %>%
  unname() %>%
  unlist()

for (directory in directories) {
  input_directory <- paste0(input_file_path, directory, "/")

  input_files <- list.files(input_directory)

  for (input_file in input_files) {
    temp_input_file_path <- paste0(input_directory, input_file)

    temp_input_file <- input_file %>%
      str_replace_all(" ", "_")

    gene_list <- fread(temp_input_file_path) %>%
      as.list() %>%
      unname() %>%
      unlist() %>%
      as.list()

    if (directory != "diseaseRelatedGenes") {
      gene_list <- gene_list[!gene_list %in% genesToRemove]
    }

    gene_list <-
      paste(directory, paste(temp_input_file, paste(gene_list, collapse = " ")), sep = "_")

    magma_input <- append(gene_list, magma_input)
  }
}

magma_input <- as.list(magma_input) %>%
  unlist() %>%
  as.list()

fwrite(
  magma_input,
  paste0(output_file_path, "geneSetsExcludingKnownDiseaseGenes.txt"),
  sep = "\n",
  quote = FALSE
)
