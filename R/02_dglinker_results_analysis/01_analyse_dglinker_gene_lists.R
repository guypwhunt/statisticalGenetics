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
  paste0("data/03_dglinker_results_analysis/01_model_results/",
         input_file_name,
         "/")

file_paths <-
  list.dirs(path = input_file_path,
            full.names = FALSE,
            recursive = FALSE)


i <- file_paths[1]

model_results <- lapply(file_paths, function(i) {
  try({
    file_path <- paste0(input_file_path,
                        i,
                        "/model_results/")

    filename <-
      list.files(path = file_path, pattern = "model_performance")
    filename <- paste0(file_path, filename)

    df <- fread(filename) %>% as.data.frame()

    colnames(df)[2] <- i

    return(df)
  })
})


model_results <- do.call("cbind", model_results)

colnames(model_results)[1] <- "metric"

model_results <-
  model_results[,!grepl("V1", colnames(model_results))]

dir.create(output_file_path,
           showWarnings = FALSE,
           recursive = TRUE)

fwrite(model_results,
       paste0(output_file_path, "model_summary.csv"),
       row.names = FALSE)
