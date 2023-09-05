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

# disease_associated_genes <-
#   c(
#     "MOBP",
#     "NEK1",
#     "TNIP1",
#     "ERGIC1",
#     "HLA-DQB1",
#     "PTPRN2",
#     "C9orf72",
#     "TBK1",
#     "KIF5A",
#     "COG3",
#     "SCFD1",
#     "SLC9A8",
#     "CFAP410",
#     "SOD1",
#     "UNC13A"
#   )

disease_associated_genes <-
  paste0(input_file_path,
         "MONDO_0004976_associations_export_als.tsv") %>%
  fread() %>%
  as.data.frame() %>%
  filter(pValue < 5 * 10 ^ -8) %>%
  # filter(
  #   traitName %in% c(
  #     "Amyotrophic lateral sclerosis",
  #     "Amyotrophic lateral sclerosis (sporadic)"
  #   )
  # ) %>%
  dplyr::select(mappedGenes) %>%
  unique() %>%
  na.omit() %>%
  as.list() %>%
  unname() %>%
  unlist() %>%
  str_split(",") %>%
  unlist() %>%
  str_split(";") %>%
  unlist() %>%
  unique() %>%
  trimws()

disease_associated_genes2 <-
  paste0(
    input_file_path,
    "gwas-association-downloaded_2023-08-02-EFO_0001357_sporadic.tsv"
  ) %>%
  fread() %>%
  as.data.frame() %>%
  filter(`P-VALUE` < 5 * 10 ^ -8) %>%
  #filter(`MAPPED_TRAIT` == "sporadic amyotrophic lateral sclerosis") %>%
  dplyr::select("REPORTED GENE(S)") %>%
  unique() %>%
  na.omit() %>%
  as.list() %>%
  unname() %>%
  unlist() %>%
  str_split(",") %>%
  unlist() %>%
  str_split(";") %>%
  unlist() %>%
  unique() %>%
  trimws()

disease_associated_genes <-
  append(disease_associated_genes, disease_associated_genes2) %>% unique()

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
