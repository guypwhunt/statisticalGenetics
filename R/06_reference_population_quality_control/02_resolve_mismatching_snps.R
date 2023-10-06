library(data.table)
library(dplyr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

input_file_name <- "als"

input_file_path <-
  paste0("data/08_quality_controlled_reference_data/",
         input_file_name,
         "/")

valid <- fread(paste0(input_file_path, "ref.valid.sample"))

dat <-
  fread(paste0(input_file_path, "ref.QC.sexcheck"))[FID %in% valid$FID]

fwrite(dat[STATUS == "OK", c("FID", "IID")], paste0(input_file_path, "ref.QC.valid"), sep =
         "\t")
