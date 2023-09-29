library(data.table)
library(dplyr)
library(stringr)

input_file_path <-
  "data/01_geneLists/08_geneListsInEntrezIdFormat/"

output_file_path <-
  "data/01_geneLists/09_magmaInput/"

dir.create(output_file_path, showWarnings = FALSE)

directories <- sub('.*\\/', '', list.dirs(input_file_path))

directories <- directories[directories != ""]

magma_input <- list()

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
  paste0(output_file_path, "geneSetsIncludingKnownDiseaseGenes.txt"),
  sep = "\n",
  quote = FALSE
)
