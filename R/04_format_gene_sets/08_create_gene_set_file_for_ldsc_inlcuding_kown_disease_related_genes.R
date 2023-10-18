library(data.table)
library(dplyr)
library(stringr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

# input_file_name <- "als"

input_file_path <-
  paste0("data/06_format_gene_sets/",
         input_file_name,
         "/01_gene_sets/")

output_file_path <-
  paste0("data/06_format_gene_sets/",
         input_file_name,
         "/03_ldscInput/geneSetsIncludingKnownDiseaseGenes/")

dir.create(output_file_path, showWarnings = FALSE)

directories <- sub('.*\\/', '', list.dirs(input_file_path))

directories <- directories[directories != ""]

for (directory in directories) {
  input_directory <- paste0(input_file_path, directory, "/")

  input_files <- list.files(input_directory)

  output_directory <- paste0(output_file_path, directory, "/")

  dir.create(output_directory, showWarnings = FALSE)

  for (input_file in input_files) {
    temp_input_file_path <- paste0(input_directory, input_file)

    temp_input_file <- input_file %>%
      str_replace_all(" ", "_")

    gene_list <- fread(temp_input_file_path) %>%
      as.list() %>%
      unname() %>%
      unlist()

    fwrite(
      as.data.frame(gene_list),
      paste0(output_directory, gsub(" ", "_", gsub("csv", "txt", input_file))),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\n"
    )

  }
}
