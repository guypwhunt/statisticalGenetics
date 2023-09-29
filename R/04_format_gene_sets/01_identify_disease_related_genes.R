library(data.table)
library(tidyverse)
library(dplyr)
library(biomaRt)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

# input_file_name <- "als"

input_file_path <-
  paste0("data/01_data_input/05_gwas_associated_genes/",
         input_file_name,
         "/")

output_file_path <-
  paste0("data/05_disease_associated_genes/", input_file_name, "/")

ensembl <-
  useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gwas_disease_associated_genes <-
  paste0(input_file_path, list.files(input_file_path)) %>%
  fread(header = FALSE) %>%
  as.list() %>%
  unlist() %>%
  unique()

gwas_disease_associated_genes <- getBM(
  attributes = c('hgnc_symbol'),
  filters = 'hgnc_symbol',
  values = gwas_disease_associated_genes,
  mart = ensembl
) %>%
  as.list() %>%
  unlist() %>%
  unique()

length(gwas_disease_associated_genes)

dir.create(output_file_path,
           showWarnings = FALSE,
           recursive = TRUE)

gwas_disease_associated_genes %>%
  as.data.frame() %>%
  fwrite(
    paste0(output_file_path, "gwasAssociatedGenes.csv"),
    row.names = FALSE,
    col.names = FALSE
  )

disease_associated_genes <-
  paste0("data/01_data_input/01_disease_related_genes/",
         input_file_name,
         ".csv") %>%
  fread(header = FALSE) %>%
  as.list() %>%
  unlist() %>%
  unname()

disease_associated_genes <- getBM(
  attributes = c('hgnc_symbol'),
  filters = 'hgnc_symbol',
  values = disease_associated_genes,
  mart = ensembl
) %>%
  as.list() %>%
  unlist() %>%
  unique()

combined_disease_gene_list <- disease_associated_genes %>%
  append(gwas_disease_associated_genes) %>%
  unique() %>%
  as.data.frame()

fwrite(
  combined_disease_gene_list,
  paste0(output_file_path, "gwasAndDiseaseAssociatedGenes.csv"),
  row.names = FALSE,
  col.names = FALSE
)

disease_associated_genes %>%
  as.data.frame() %>%
  fwrite(
    paste0(output_file_path, "diseaseAssociatedGenes.csv"),
    row.names = FALSE,
    col.names = FALSE
  )
