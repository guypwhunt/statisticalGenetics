library(data.table)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(ggvenn)

input_file_path <-
  "data/01_geneLists/07_diseaseGwasResults/gwasCatalogue/"

output_file_path <-
  "data/01_geneLists/07_diseaseGwasResults/genesToRemove/"

ensembl <-
  useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

disease_associated_genes <-
  c(
    "MOBP",
    "NEK1",
    "TNIP1",
    "ERGIC1",
    "HLA-DQB1",
    "PTPRN2",
    "C9orf72",
    "TBK1",
    "KIF5A",
    "COG3",
    "SCFD1",
    "SLC9A8",
    "CFAP410",
    "SOD1",
    "UNC13A"
  )

disease_associated_genes <- getBM(
  attributes = c('hgnc_symbol'),
  filters = 'hgnc_symbol',
  values = disease_associated_genes,
  mart = ensembl
) %>%
  as.list() %>%
  unlist() %>%
  unique()

length(disease_associated_genes)

dir.create(output_file_path,
           showWarnings = FALSE)

write.csv(
  disease_associated_genes,
  paste0(output_file_path, "gwasCatalogueGenes.csv"),
  row.names = FALSE
)

original_disease_gene_list <-
  read.csv("data/01_geneLists/01_knownDiseaseRelatedGenes/combinedGeneList.csv") %>%
  as.list() %>%
  unlist() %>%
  unname()

combined_disease_gene_list <- original_disease_gene_list %>%
  append(disease_associated_genes) %>%
  unique()

write.csv(
  combined_disease_gene_list,
  paste0(output_file_path, "gwasCataloguePlusOrginalGenes.csv"),
  row.names = FALSE
)

combined_gene_list <- list(original = original_disease_gene_list,
                           gwas = disease_associated_genes)

ggvenn(combined_gene_list,
       set_name_size = 12,
       text_size = 12)
