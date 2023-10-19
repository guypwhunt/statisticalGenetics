library(data.table)
library(dplyr)
library(tools)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]
mfa <- commandArgs(trailingOnly = TRUE)[2] %>% as.numeric() %>%
  format(scientific = FALSE)

# input_file_name <- "als"
# mfa <- 0.0025 %>%
#   format(scientific = FALSE)

input_file_path <-
  paste0("data/12_ldsc_results/", input_file_name, "/",
         mfa,
         "/")

output_file_path <-
  paste0("data/13_consolidated_ldsc_results/",
         input_file_name,
         "/",
         mfa,
         "/")

dir.create(output_file_path, showWarnings = FALSE, recursive = TRUE)

parentDirectories <-
  list.dirs(input_file_path, recursive = FALSE, full.names = FALSE)

for (parentDirectory in parentDirectories) {
  input_parent_file_path <-
    paste0(input_file_path, parentDirectory, "/")

  childDirectories <-
    list.dirs(input_parent_file_path,
              recursive = FALSE,
              full.names = FALSE)

  for (childDirectory in childDirectories) {
    input_child_file_path <-
      paste0(input_parent_file_path, childDirectory, "/")

    file_names <-
      list.files(input_child_file_path, pattern = ".results")

    for (file_name in file_names) {
      ldsc_result <-
        paste0(input_child_file_path, file_name) %>%
        fread()

      ldsc_result$gene_set <-
        paste0(childDirectory, "_", file_path_sans_ext(file_name))

      if (exists("consolidated_results")) {
        consolidated_results <-
          rbind(consolidated_results, ldsc_result)
      } else {
        consolidated_results <- ldsc_result
      }
    }
  }

  consolidated_results %>%
    fwrite(paste0(output_file_path, parentDirectory, ".csv"))

  rm("consolidated_results")
}
