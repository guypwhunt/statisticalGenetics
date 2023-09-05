library(data.table)
library(biomaRt)
library(dplyr)

input_file_path <-
  "data/01_geneLists/07_diseaseGwasResults/genesToRemove/"

original_genes_input_file_path <-
  "data/01_geneLists/01_knownDiseaseRelatedGenes/combinedGeneList.csv"

output_file_path <-
  "data/01_geneLists/08_geneListsInEntrezIdFormat/diseaseRelatedGenes/"

dir.create(output_file_path,
           showWarnings = FALSE,
           recursive = TRUE)

input_file_paths <-
  c(original_genes_input_file_path,
    paste0(input_file_path, list.files(input_file_path)))

ensembl <-
  useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

for (file_path in input_file_paths) {
  disease_related_genes <-
    fread(file_path) %>%
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

  fwrite(disease_related_genes,
         paste0(output_file_path,
                sub('.*\\/', '', file_path)),
         row.names = FALSE)
}
