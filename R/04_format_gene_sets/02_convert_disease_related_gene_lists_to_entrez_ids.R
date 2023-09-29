library(data.table)
library(biomaRt)
library(dplyr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

# input_file_name <- "als"

input_file_path <-
  paste0("data/05_disease_associated_genes/", input_file_name, "/")

output_file_path <-
  paste0("data/06_format_gene_sets/",
         input_file_name,
         "/01_gene_sets/diseaseAssociatedGenes/")

dir.create(output_file_path,
           showWarnings = FALSE,
           recursive = TRUE)

input_file_paths <-
  paste0(input_file_path, list.files(input_file_path))

ensembl <-
  useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

for (file_path in input_file_paths) {
  disease_related_genes <-
    fread(file_path, header = FALSE) %>%
    as.list() %>%
    unname() %>%
    unlist()

  x <- listAttributes(ensembl)

  disease_related_genes <- getBM(
    attributes = c('entrezgene_id'),
    filters = 'hgnc_symbol',
    values = disease_related_genes,
    mart = ensembl
  ) %>%
    unique()

  fwrite(
    disease_related_genes,
    paste0(output_file_path,
           sub('.*\\/', '', file_path)),
    row.names = FALSE,
    col.names = FALSE
  )
}
