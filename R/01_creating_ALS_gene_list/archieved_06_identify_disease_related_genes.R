library(data.table)
library(tidyverse)
library(dplyr)
library(R.utils)
library(biomaRt)
library(snow)
library(parallel)

disease_associated_genes <-
  fread(
    "data/geneLists/06_diseaseGwasResults/summaryStats/GCST90013429_buildGRCh37.tsv"
  ) %>%
  as.data.frame() %>%
  arrange(p_value) %>%
  filter(p_value < 0.05) %>%
  arrange(chromosome, base_pair_location)

#####################################
mart <- useEnsembl(biomart = "snps",
                   dataset = "hsapiens_snp")

attributes <- listAttributes(mart)

filters <-
  listFilters(mart, c("name", "description", "fullDescription"))

chromosome_location <-
  data.frame(CHR = disease_associated_genes[, 1],
             START = disease_associated_genes[, 2],
             END = disease_associated_genes[, 2]) %>%
  apply(1, paste, collapse = ":") %>% unique()

cl <- makeCluster(detectCores() - 1, type = "SOCK")

clusterExport(cl, c("mart", "getBM"))

ensemblGeneIds <-
  clusterApply(cl, chromosome_location, function(i) {
    try({
      getBM(
        attributes = c('ensembl_gene_stable_id'),
        filters = c('chromosomal_region'),
        values = i,
        mart = mart
      )
    })
  }) %>%
  unlist() %>%
  unique() %>%
  na.omit()

stopCluster(cl)

ensemblMart <-
  useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

geneSymbols <-
  getBM(
    filters = "ensembl_gene_id",
    attributes = c('hgnc_symbol'),
    values = ensemblGeneIds,
    mart = ensemblMart
  ) %>%
  unlist() %>%
  unique() %>%
  na.omit()

print(geneSymbols)

backupGeneSymbols <- geneSymbols
