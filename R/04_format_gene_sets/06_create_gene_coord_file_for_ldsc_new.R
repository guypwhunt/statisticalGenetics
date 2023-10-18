library(data.table)
library(dplyr)
library(stringr)
library(biomaRt)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

# input_file_name <- "als"

chromosomes <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
                 "MT", "X", "Y")

ensembl <-
  useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

input_file_path <-
  paste0("data/06_format_gene_sets/",
         input_file_name,
         "/01_gene_sets/")

output_file_path <-
  paste0("data/06_format_gene_sets/",
         input_file_name,
         "/03_ldscInput/")

dir.create(output_file_path, showWarnings = FALSE)

directories <- sub('.*\\/', '', list.dirs(input_file_path))

directories <- directories[directories != ""]

ldsc_input <- list()

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
      unlist()

    ldsc_input <- append(ldsc_input, gene_list)
  }
}

ldsc_input <- ldsc_input %>% unlist() %>% unique()

attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)

gene_coord_file <- getBM(
  attributes = c(
    'entrezgene_id',
    'chromosome_name',
    'start_position',
    'end_position'
  ),
  filters = 'entrezgene_id',
  values = ldsc_input,
  mart = ensembl
) %>%
  unique() %>%
  rename(GENE = entrezgene_id,
         CHR = chromosome_name,
         START = start_position,
         END = end_position) %>%
  arrange(GENE, CHR) %>%
  filter(CHR %in% chromosomes)

gene_coord_file <- gene_coord_file[!duplicated(gene_coord_file$GENE),]

fwrite(
  gene_coord_file,
  paste0(output_file_path, "geneCoordFile.txt"),
  sep = "\t",
  quote = FALSE
)
