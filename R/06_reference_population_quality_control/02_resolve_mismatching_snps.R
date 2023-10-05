library(data.table)
library(dplyr)
library(stringr)
library(magrittr)

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

bim <- fread("EUR.bim") %>%
  setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
  .[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]

height <- fread("Height.QC.gz") %>%
  .[,c("A1","A2"):=list(toupper(A1), toupper(A2))]

qc <- fread("EUR.QC.snplist", header=F)

# Merge summary statistic with target
info <- merge(bim, height, by=c("SNP", "CHR", "BP")) %>%
  # And filter out QCed SNPs
  .[SNP %in% qc[,V1]]

# Function for calculating the complementary allele
complement <- function(x){
  switch (x,
          "A" = "T",
          "C" = "G",
          "T" = "A",
          "G" = "C",
          return(NA)
  )
}
# Get SNPs that have the same alleles across base and target
info.match <- info[A1 == B.A1 & A2 == B.A2, SNP]
# Identify SNPs that are complementary between base and target
com.snps <- info[sapply(B.A1, complement) == A1 &
                   sapply(B.A2, complement) == A2, SNP]
# Now update the bim file
bim[SNP %in% com.snps, c("B.A1", "B.A2") :=
      list(sapply(B.A1, complement),
           sapply(B.A2, complement))]

# identify SNPs that need recoding
recode.snps <- info[B.A1==A2 & B.A2==A1, SNP]
# Update the bim file
bim[SNP %in% recode.snps, c("B.A1", "B.A2") :=
      list(B.A2, B.A1)]

# identify SNPs that need recoding & complement
com.recode <- info[sapply(B.A1, complement) == A2 &
                     sapply(B.A2, complement) == A1, SNP]
# Now update the bim file
bim[SNP %in% com.recode, c("B.A1", "B.A2") :=
      list(sapply(B.A2, complement),
           sapply(B.A1, complement))]
# Write the updated bim file
fwrite(bim[,c("SNP", "B.A1")], "EUR.a1", col.names=F, sep="\t")

mismatch <- bim[!(SNP %in% info.match |
                    SNP %in% com.snps |
                    SNP %in% recode.snps |
                    SNP %in% com.recode), SNP]
write.table(mismatch, "EUR.mismatch", quote=F, row.names=F, col.names=F)
