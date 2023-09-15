library(data.table)
library(tidyverse)
library(dplyr)
library(stringr)
library(readxl)
library(biomaRt)
library(ggvenn)

input_file_path <-
  "data/01_data_input/raw_inputs/"

output_file_path <-
  "data/01_data_input/01_disease_related_genes/"

ensembl <-
  useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

disgenet_gene_list <-
  fread(paste0(input_file_path, "Disgenet.tsv")) %>%
  as.data.frame() %>%
  filter(Score_gda >= 0.3) %>%
  dplyr::select(Gene) %>%
  as.list() %>%
  unlist() %>%
  str_trim() %>%
  unique()

disgenet_gene_list <- getBM(
  attributes = c('hgnc_symbol'),
  filters = 'hgnc_symbol',
  values = disgenet_gene_list,
  mart = ensembl
) %>%
  as.list() %>%
  unlist() %>%
  unique()

head(disgenet_gene_list)
length(disgenet_gene_list)

clinvar_clinical_significance_column_name <-
  "Clinical significance (Last reviewed)"

clinvar_gene_list <-
  fread(paste0(input_file_path, "clinvar.txt")) %>%
  as.data.frame() %>%
  filter(grepl(
    "Pathogenic",
    `Clinical significance (Last reviewed)`,
    ignore.case = FALSE
  ))  %>%
  filter(
    grepl("Amyotrophic lateral sclerosis",
          `Condition(s)`,
          ignore.case = TRUE) | grepl("motor neurone disease",
                                      `Condition(s)`,
                                      ignore.case = TRUE) |
      grepl("motor neuron disease",
            `Condition(s)`,
            ignore.case = TRUE)
  ) %>%
  dplyr::select("Gene(s)") %>%
  as.list() %>%
  unname() %>%
  unlist() %>%
  str_split("\\|") %>%
  unlist() %>%
  str_trim() %>%
  unique()

clinvar_gene_list <- getBM(
  attributes = c('hgnc_symbol'),
  filters = 'hgnc_symbol',
  values = clinvar_gene_list,
  mart = ensembl
) %>%
  as.list() %>%
  unlist() %>%
  unique()

head(clinvar_gene_list)
length(clinvar_gene_list)

alsod_gene_list <-
  read_excel(paste0(input_file_path, "alsod.xlsx")) %>%
  as.data.frame() %>%
  filter(Category != 'Tenuous' |
           Category == 'Unassigned') %>%
  dplyr::select("Gene symbol") %>%
  as.list() %>%
  unname() %>%
  unlist() %>%
  str_trim() %>%
  unique()

alsod_gene_list <- getBM(
  attributes = c('hgnc_symbol'),
  filters = 'hgnc_symbol',
  values = alsod_gene_list,
  mart = ensembl
) %>%
  as.list() %>%
  unlist() %>%
  unique()

head(alsod_gene_list)
length(alsod_gene_list)

omim_gene_list <-
  read_excel(paste0(input_file_path, "OMIM.xlsx"), skip = 4) %>%
  as.data.frame() %>%
  filter(
    grepl("Amyotrophic lateral sclerosis",
          `Phenotype`,
          ignore.case = TRUE) | grepl("motor neurone disease",
                                      `Phenotype`,
                                      ignore.case = TRUE) |
      grepl("motor neuron disease",
            `Phenotype`,
            ignore.case = TRUE)
  ) %>%
  dplyr::select("Gene/Locus") %>%
  as.list() %>%
  unname() %>%
  unlist() %>%
  str_split("\\,") %>%
  unlist() %>%
  str_trim() %>%
  unique()

omim_gene_list <- getBM(
  attributes = c('hgnc_symbol'),
  filters = 'hgnc_symbol',
  values = omim_gene_list,
  mart = ensembl
) %>%
  as.list() %>%
  unlist() %>%
  unique()

head(omim_gene_list)
length(omim_gene_list)

number_of_genes_per_list <- data.frame(list = as.character(),
                                       numberOfGenes = as.integer())
combined_gene_list <- list(
  omim = omim_gene_list,
  alsod = alsod_gene_list,
  clinvar = clinvar_gene_list,
  disgenet = disgenet_gene_list
)

ggvenn(combined_gene_list)
ggsave(
  paste0(input_file_path,
         "summaryStats/vennDiagram.png"),
  dpi = 600,
  width = 12,
  height = 8,
  units = c("cm")
)


for (list_number in seq(length(combined_gene_list))) {
  number_of_genes_per_list[list_number, ] <-
    c(names(combined_gene_list)[list_number],
      length(combined_gene_list[[list_number]]))
}

dir.create(paste0(input_file_path, "summaryStats"), showWarnings = FALSE)

write.csv(
  number_of_genes_per_list,
  paste0(input_file_path, "summaryStats/numberOfGenesPerList.csv")
)

number_of_common_genes_per_list <-
  data.frame(
    listOne = as.character(),
    listTwo = as.character(),
    numberOfCommonGenes = as.integer()
  )

row_number <- 1

for (list_number_one in seq(length(combined_gene_list))) {
  for (list_number_two in seq(length(combined_gene_list))) {
    number_of_common_genes_per_list[row_number,] <-
      c(names(combined_gene_list)[list_number_one],
        names(combined_gene_list)[list_number_two],
        length(intersect(
          combined_gene_list[[list_number_one]], combined_gene_list[[list_number_two]]
        )))

    row_number <- row_number + 1
  }
}

write.csv(
  number_of_common_genes_per_list,
  paste0(
    input_file_path,
    "summaryStats/numberOfCommonGenesPerList.csv"
  )
)

combined_gene_list <- combined_gene_list %>%
  unlist() %>%
  unname() %>%
  unique() %>%
  str_replace_all("@", "")

write.csv(combined_gene_list,
          paste0(output_file_path, "als.csv"),
          row.names = FALSE)
