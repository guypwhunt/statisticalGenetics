library(GO.db)
library(dplyr)
library(tidyverse)

file_path <- "data/01_geneLists/02_geneEnrichmentResults"

input_file_path <-
  "data/01_geneLists/02_geneEnrichmentResults/"

output_file_path <-
  "data/01_geneLists/03_significantGeneEnrichmentResults/"

go_hierarchy <- as.data.frame(GO.db::GOBPANCESTOR)
go_hierarchy <-
  rbind(go_hierarchy, as.data.frame(GO.db::GOMFANCESTOR))
go_hierarchy <-
  rbind(go_hierarchy, as.data.frame(GO.db::GOCCANCESTOR))

# Uncomment to find parents
# go_hierarchy <- go_hierarchy[, c(2,1)]


go_file_names <-
  list.files(
    path = input_file_path,
    pattern = "GO_",
    all.files = FALSE,
    full.names = FALSE
  )

go_file_names <- go_file_names[grep("2023", go_file_names)]

go_file_paths <- paste0(input_file_path, go_file_names)

go_data_frames <- lapply(go_file_paths, function(i) {
  read.csv(i, header = TRUE)
})

dir.create(
  output_file_path,
  showWarnings = FALSE,
  recursive = TRUE
)

for (go_data_frame_number in seq(length(go_data_frames))) {
  temp_go_data_frame <-
    go_data_frames[[go_data_frame_number]] %>%
    dplyr::select(-c(X, Old.P.value, Old.Adjusted.P.value)) %>%
    filter(Adjusted.P.value < 0.05) %>%
    mutate(ID = str_remove(str_split(Term, "\\(", simplify = T)[, 2], "\\)"))

  message(nrow(temp_go_data_frame))

  temp_go_hierarchy <-
    go_hierarchy[go_hierarchy[, 1] %in% temp_go_data_frame$ID,]

  message(nrow(temp_go_data_frame))

  go_data_frames[[go_data_frame_number]] <- temp_go_data_frame

  write.csv(
    temp_go_data_frame,
    file = paste0(
      output_file_path,
      go_file_names[go_data_frame_number]
    ),
    row.names = FALSE
  )
}
