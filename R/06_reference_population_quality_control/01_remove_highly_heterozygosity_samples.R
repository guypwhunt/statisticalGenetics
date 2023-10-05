library(data.table)
library(dplyr)
library(stringr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

input_file_name <- "als"

input_file_path <-
  paste0("data/08_quality_controlled_reference_data/",
         input_file_name,
         "/")

output_file_path <-
  paste0("data/08_quality_controlled_reference_data/",
         input_file_name,
         "/")

ref_data <- fread(paste0(input_file_path, "ref.QC.het"))

valid <- ref_data[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)]

fwrite(valid[,c("FID","IID")], paste0(output_file_path,"ref.valid.sample"), sep="\t")
