library(data.table)
library(dplyr)
library(stringr)
library(magrittr)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]
mfa <- commandArgs(trailingOnly = TRUE)[2] %>% as.numeric() %>%
  format(scientific = FALSE)

# input_file_name <- "als"

reference_population_input_file_path <-
  paste0("data/01_data_input/07_reference_population/",
         input_file_name,
         "/")

gwas_input_file_path <-
  paste0(
    "data/07_quality_controlled_gwas_data/",
    input_file_name,
    "/",
    mfa,
    "/",
    "/quality_controlled_gwas_data.tsv"
  )

qcd_snps_file_path <-
  paste0("data/08_quality_controlled_reference_data/",
         input_file_name,
         "/",
         mfa,
         "/ref.QC.snplist")

output_file_path <-
  paste0("data/08_quality_controlled_reference_data/",
         input_file_name,
         "/",
         mfa,
         "/")

bim_data <-
  paste0(
    reference_population_input_file_path,
    list.files(reference_population_input_file_path, pattern = ".bim")
  ) %>%
  fread() %>%
  setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
  .[, c("B.A1", "B.A2") := list(toupper(B.A1), toupper(B.A2))]

gwas_data <- fread(gwas_input_file_path) %>%
  .[, c("effect_allele", "other_allele") := list(toupper(effect_allele), toupper(other_allele))]

qc <- fread(qcd_snps_file_path, header = FALSE) %>% c() %>%
  unname() %>%
  unlist()


info <-
  merge(bim_data,
        gwas_data,
        by.x = c("SNP", "CHR", "BP"),
        by.y = c("rsid", "chromosome", "base_pair_location")) %>%
  .[SNP %in% qc]

complement <- function(x) {
  switch (
    x,
    "A" = "T",
    "C" = "G",
    "T" = "A",
    "G" = "C",
    return(NA)
  )
}


info.match <- info[effect_allele == B.A1 & other_allele == B.A2, SNP]

com.snps <- info[sapply(B.A1, complement) == effect_allele &
                   sapply(B.A2, complement) == other_allele, SNP]

bim_data <- bim_data[SNP %in% com.snps, c("B.A1", "B.A2") :=
                       list(sapply(B.A1, complement),
                            sapply(B.A2, complement))]

recode.snps <- info[B.A1 == other_allele & B.A2 == effect_allele, SNP]

bim_data <- bim_data[SNP %in% recode.snps, c("B.A1", "B.A2") :=
           list(B.A2, B.A1)]

com.recode <- info[sapply(B.A1, complement) == other_allele &
                     sapply(B.A2, complement) == effect_allele, SNP]

bim_data <- bim_data[SNP %in% com.recode, c("B.A1", "B.A2") :=
           list(sapply(B.A2, complement),
                sapply(B.A1, complement))]

fwrite(bim_data[, c("SNP", "B.A1")], paste0(output_file_path,"ref.a1"), col.names = FALSE, sep = "\t")

mismatch <- bim_data[!(SNP %in% info.match |
                         SNP %in% com.snps |
                         SNP %in% recode.snps |
                         SNP %in% com.recode), SNP]
fwrite(
  x = as.data.frame(mismatch),
  file = paste0(output_file_path, "ref.mismatch"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
