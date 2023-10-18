library(data.table)
library(dplyr)
library(stringr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

mfa <- commandArgs(trailingOnly = TRUE)[2] %>% as.numeric() %>%
  format(scientific = FALSE)

# input_file_name <- "als"

input_file_path <-
  paste0("data/01_data_input/06_gwas_summary_statistics/",
         input_file_name,
         "/")

output_file_path <-
  paste0("data/07_quality_controlled_gwas_data/",
         input_file_name,
         "/",
         mfa,
         "/")

dir.create(output_file_path,
           showWarnings = FALSE,
           recursive = TRUE)

input_file_path <-
  paste0(input_file_path, list.files(input_file_path))

gwas_data <- fread(input_file_path) %>%
  filter(effect_allele_frequency > mfa)

totalRows <- nrow(gwas_data)

if(length(unique(gwas_data$rsid)) != totalRows) {
  rsids <- c()

  for(row in seq(totalRows)) {
    temp_gwas_data <- gwas_data[row , ]

    if (!temp_gwas_data$rsid %in% rsids) {
      rsids <- append(rsids, temp_gwas_data$rsid)

      if (exists("filtered_gwas_data")) {
        filtered_gwas_data <- rbind(filtered_gwas_data, temp_gwas_data)
      } else {
        filtered_gwas_data <- temp_gwas_data
      }
    }
  }
  gwas_data <- filtered_gwas_data
}

if ("info" %in% colnames(gwas_data)) {
  gwas_data <- gwas_data %>%
    filter(info > 0.8)
}

filter_condition <-
  !(
    (gwas_data$effect_allele == "A" &
       gwas_data$other_allele == "T") |
      (gwas_data$effect_allele == "T" &
         gwas_data$other_allele == "A") |
      (gwas_data$effect_allele == "G" &
         gwas_data$other_allele == "C") |
      (gwas_data$effect_allele == "C" &
         gwas_data$other_allele == "G")
  )

gwas_data <- subset(gwas_data, filter_condition)


fwrite(gwas_data,
       paste0(output_file_path, "quality_controlled_gwas_data.tsv"),
       sep = "\t")
