library(data.table)
library(tidyverse)
library(dplyr)
library(ggvenn)
library(stringr)

input_file_path <-
  "data/01_geneLists/06_dglinker_results"

file_paths <-
  list.dirs(path = input_file_path,
            full.names = TRUE,
            recursive = FALSE) %>%
  paste0("/model_results/")

model_results <- lapply(file_paths, function(i) {
  try({
    filename <- list.files(path = i, pattern = "results")
    filename <- paste0(i, filename)

    df <- fread(filename) %>% as.data.frame() %>%
      dplyr::filter(`Association type` != "Not-associated") %>%
      dplyr::select('Gene name') %>% as.list() %>% unlist() %>% as.list() %>%
      unique() %>% unlist()

    return(df)
  })
})

names(model_results) <- sapply(str_split(file_paths, "/"), `[[`, 4)

model_results <- model_results[c(1,2,4,5)]

ggvenn(model_results, set_name_size = 2)
