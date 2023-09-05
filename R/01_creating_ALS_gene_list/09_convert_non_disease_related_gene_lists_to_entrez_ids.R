library(data.table)
library(biomaRt)
library(dplyr)

input_file_path <-
  "data/01_geneLists/05_geneEnrichmentTermGenes/"

output_file_path <-
  "data/01_geneLists/08_geneListsInEntrezIdFormat/"

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

    gene_list <- fread(temp_input_file_path) %>%
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

    fwrite(gene_list,
           paste0(output_directory,
                  input_file),
           row.names = FALSE)

  }
}

dglinker_genes_input_file_path <-
  "data/01_geneLists/06_dglinker_results/allPublicationDataBases/model_results/"

dglinker_file_name <- list.files(dglinker_genes_input_file_path,
                                 pattern = "results")

dglinker_file_path <- paste0(dglinker_genes_input_file_path,
                             dglinker_file_name)

dglinker_data <- fread(dglinker_file_path) %>%
  as.data.frame() %>%
  filter(`Association type` != "Not-associated") %>%
  dplyr::select(`Gene name`) %>% as.list() %>% unname() %>% unlist()

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

fwrite(dglinker_data,
       paste0(output_directory,
              "dglinker.csv"),
       row.names = FALSE)
