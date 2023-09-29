library(data.table)
library(tidyverse)
library(dplyr)
library(ggvenn)
library(stringr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

# input_file_name <- "als"

input_file_path <-
  paste0("data/01_data_input/03_dglinker_results/",
         input_file_name,
         "/")

output_file_path <-
  paste0("data/03_dglinker_results_analysis/02_dglinker_gene_list/",
         input_file_name,
         "/")

model_summary_file_path <-
  paste0("data/03_dglinker_results_analysis/01_model_results/",
         input_file_name,
         "/")

model_results <-
  fread(paste0(model_summary_file_path, "model_summary.csv")) %>%
  as.data.frame() %>%
  column_to_rownames(var = "metric")

best_model <-
  colnames(model_results)[max.col(model_results[row.names(model_results) == "AUC of the Model",])]

input_file_path <-
  paste0(input_file_path, best_model, "/model_results/")

input_file_path <-
  paste0(input_file_path,
         list.files(input_file_path, pattern = "results"))

df <- fread(input_file_path) %>% as.data.frame() %>%
  dplyr::filter(`Association type` != "Not-associated") %>%
  dplyr::select('Gene name') %>% as.list() %>% unlist() %>% as.list() %>%
  unique() %>% unlist() %>%
  as.data.frame()

dir.create(output_file_path,
           showWarnings = FALSE,
           recursive = TRUE)

fwrite(
  df,
  paste0(output_file_path, "dglinker_gene_list.csv"),
  row.names = FALSE,
  col.names = FALSE
)
