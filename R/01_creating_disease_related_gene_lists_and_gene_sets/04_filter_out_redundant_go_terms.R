library(dplyr)
library(tidyverse)
library(rrvgo)

input_file_path <-
  "data/01_geneLists/03_significantGeneEnrichmentResults/"

output_file_path <-
  "data/01_geneLists/04_significantNonRedundantGeneEnrichmentResults/"

go_file_names <-
  list.files(
    path = input_file_path,
    pattern = "GO",
    all.files = FALSE,
    full.names = FALSE
  )

go_file_name_to_summarise <-
  go_file_names[grep("Biological", go_file_names)]

go_file_name_not_to_summarise <-
  go_file_names[go_file_names != go_file_name_to_summarise]

go_path_name_to_summarise <-
  paste0(input_file_path, "/", go_file_name_to_summarise)

go_data_frame <-
  read.csv(go_path_name_to_summarise, header = TRUE)

dir.create(output_file_path,
           showWarnings = FALSE,
           recursive = TRUE)

simMatrix <- calculateSimMatrix(go_data_frame$ID,
                                orgdb = "org.Hs.eg.db",
                                ont = "BP",
                                method = "Rel")


scores <-
  setNames(-log10(go_data_frame$Adjusted.P.value), go_data_frame$ID)

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold = 0.7,
                                orgdb = "org.Hs.eg.db")

heatmapPlot(
  simMatrix,
  reducedTerms,
  annotateParent = TRUE,
  annotationLabel = "parentTerm",
  fontsize = 6
)

scatterPlot(simMatrix, reducedTerms)

treemapPlot(reducedTerms)

wordcloudPlot(reducedTerms, min.freq = 1, colors = "black")

reducedGoTerms <-
  data.frame(matrix(nrow = 0, ncol = length(colnames(reducedTerms))))
colnames(reducedGoTerms) <- colnames(reducedTerms)

for (cluster in seq(unique(reducedTerms$cluster))) {
  temp_reducedTerms <- reducedTerms %>%
    filter(cluster == !!cluster) %>%
    filter(score == max(score))

  reducedGoTerms <- rbind(reducedGoTerms, temp_reducedTerms)
}

go_data_frame <-
  go_data_frame[go_data_frame$ID %in% reducedGoTerms$go,]

write.csv(go_data_frame,
          paste0(output_file_path, go_file_name_to_summarise),
          row.names = FALSE)

for(file_name in go_file_name_not_to_summarise){
  temp_df <- read.csv(paste0(input_file_path, file_name))

  write.csv(temp_df, paste0(output_file_path, file_name), row.names = FALSE)
}
