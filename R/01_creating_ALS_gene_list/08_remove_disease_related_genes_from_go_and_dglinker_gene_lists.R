library(data.table)
library(tidyverse)
library(dplyr)

genes_to_remove_input_file_path <-
  "data/01_geneLists/07_diseaseGwasResults/genesToRemove/"

go_genes_input_file_path <-
  "data/01_geneLists/05_geneEnrichmentTermGenes/"

dglinker_genes_input_file_path <-
  "data/01_geneLists/06_dglinker_results/allPublicationDataBases/model_results/"

output_file_path <-
  "data/01_geneLists/08_geneListsWithKnownDiseaseAssociatedGenesRemoved/"

dir.create(output_file_path, showWarnings = FALSE)

genes_to_remove_file_name <-
  list.files(genes_to_remove_input_file_path,
             pattern = "gwasCataloguePlusOrginalGenes")

genes_to_remove_file_path <-
  paste0(genes_to_remove_input_file_path, genes_to_remove_file_name)

genes_to_remove_data <-
  read.csv(genes_to_remove_file_path) %>% as.list() %>% unname() %>% unlist()

go_directory_names <-
  list.dirs(go_genes_input_file_path, full.names = FALSE) %>%
  na.omit()

go_directory_names <- go_directory_names[go_directory_names != ""]

for (go_directory_name in go_directory_names) {
  temp_go_genes_input_file_path <-
    paste0(go_genes_input_file_path, go_directory_name, "/")

  temp_output_file_path <-
    paste0(output_file_path, go_directory_name, "/")

  dir.create(temp_output_file_path, showWarnings = FALSE)

  temp_files <- list.files(temp_go_genes_input_file_path)

  for (temp_file in temp_files) {
    temp_data <-
      read.csv(paste0(temp_go_genes_input_file_path, temp_file)) %>%
      as.list() %>% unname() %>% unlist()

    temp_data <- temp_data[!temp_data %in% genes_to_remove_data]

    write.csv(temp_data,
              paste0(temp_output_file_path, temp_file),
              row.names = FALSE)
  }
}

dglinker_file_name <- list.files(dglinker_genes_input_file_path,
                                 pattern = "results")

dglinker_file_path <- paste0(dglinker_genes_input_file_path,
                             dglinker_file_name)

dglinker_data <- fread(dglinker_file_path) %>%
  as.data.frame() %>%
  filter(`Association type` == "Predicted") %>%
  dplyr::select(`Gene name`) %>% as.list() %>% unname() %>% unlist()

dglinker_data <-
  dglinker_data[!dglinker_data %in% genes_to_remove_data]

dglinker_output_file_path <- paste0(output_file_path, "dglinker/")

dir.create(dglinker_output_file_path, showWarnings = FALSE)

write.csv(dglinker_data,
          paste0(dglinker_output_file_path, "dglinker.csv"),
          row.names = FALSE)
