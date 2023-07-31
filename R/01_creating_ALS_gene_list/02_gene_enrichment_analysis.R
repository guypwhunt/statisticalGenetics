library(data.table)
library(tidyverse)
library(enrichR)
library(dplyr)

combined_gene_list <-
  fread("data/geneLists/combinedGeneList.csv") %>%
  as.data.frame() %>%
  dplyr::select(x) %>%
  as.list() %>%
  unlist() %>%
  unname() %>%
  c()

dbs <- listEnrichrDbs()

head(dbs)

enriched <- enrichr(combined_gene_list, dbs$libraryName)

head(enriched)

dir.create(
  "data/geneLists/enrichmentResults",
  recursive = TRUE,
  showWarnings = FALSE
)

for (n in seq(length(enriched))) {
  write.csv(enriched[[n]],
            paste0(
              "data/geneLists/enrichmentResults/",
              names(enriched)[n],
              ".csv"
            ))

  print(plotEnrich(
    enriched[[n]],
    showTerms = 20,
    numChar = 40,
    y = "Count",
    orderBy = "P.value",
    title = names(enriched)[n]
  ))
}
