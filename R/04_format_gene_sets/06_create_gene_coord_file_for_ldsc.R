library(data.table)
library(dplyr)
library(stringr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

# input_file_name <- "als"

input_file_path <-
  paste0("data/01_data_input/08_gene_locations/",
         input_file_name,
         "/")

input_file_path <-
  paste0(input_file_path, list.files(input_file_path))

output_file_path <-
  paste0("data/06_format_gene_sets/",
         input_file_name,
         "/03_ldscInput/")


fread(input_file_path, header = FALSE) %>%
  rename(
    GENE = V1,
    CHR = V2,
    START = V3,
    END = V4
  ) %>%
  select(c(GENE, CHR, START, END)) %>%
  mutate(CHR = paste0("chr", CHR)) %>%
  fwrite(paste0(output_file_path, "geneCoordFile.txt"),
         sep = "\t",
         quote = FALSE)
