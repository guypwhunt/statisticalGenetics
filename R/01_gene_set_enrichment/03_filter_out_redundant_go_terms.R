library(dplyr)
library(tidyverse)
library(rrvgo)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]

# input_file_name <- "als"

input_file_path <-
  paste0(
    "data/02_gene_set_enrichment/02_significant_gene_enrichment_results/",
    input_file_name,
    "/"
  )

output_file_path <-
  paste0(
    "data/02_gene_set_enrichment/03_significant_non_redundant_gene_enrichmentResults/",
    input_file_name,
    "/"
  )

go_file_names <-
  list.files(
    path = input_file_path,
    pattern = "GO",
    all.files = FALSE,
    full.names = FALSE
  )

#######################
dir.create(output_file_path,
           showWarnings = FALSE,
           recursive = TRUE)

for (go_file_name in go_file_names) {
  go_path_name_to_summarise <-
    paste0(input_file_path, "/", go_file_name)

  go_data_frame <-
    read.csv(go_path_name_to_summarise, header = TRUE)


  simMatrix <- calculateSimMatrix(
    go_data_frame$ID,
    orgdb = "org.Hs.eg.db",
    ont = if (go_file_name == "GO_Biological_Process_2023.csv") {
      "BP"
    } else if (go_file_name == "GO_Cellular_Component_2023.csv") {
      "CC"
    } else {
      "MF"
    },
    method = "Rel"
  )


  scores <-
    setNames(-log10(go_data_frame$Adjusted.P.value), go_data_frame$ID)

  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold = 0.7,
                                  orgdb = "org.Hs.eg.db")

  reducedGoTerms <-
    data.frame(matrix(nrow = 0, ncol = length(colnames(reducedTerms))))
  colnames(reducedGoTerms) <- colnames(reducedTerms)

  for (cluster in seq(unique(reducedTerms$cluster))) {
    temp_reducedTerms <- reducedTerms %>%
      filter(cluster == !!cluster) %>%
      arrange(desc(score)) %>%
      slice(1:3)

    reducedGoTerms <- rbind(reducedGoTerms, temp_reducedTerms)
  }

  go_data_frame <-
    go_data_frame[go_data_frame$ID %in% reducedGoTerms$go, ]

  write.csv(go_data_frame,
            paste0(output_file_path, go_file_name),
            row.names = FALSE)
}
