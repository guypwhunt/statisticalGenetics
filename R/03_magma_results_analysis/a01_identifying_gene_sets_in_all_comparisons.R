library(data.table)
library(dplyr)

input_file_path <-
  "data/02_magma/08_geneSetAnalysis/"

fileNames <- c(0, 1, 2.5, 5, 7.5, 10)

gene_set_results <- sapply(fileNames, function(fileName) {
  paste0(input_file_path, fileName, ".gsa.out") %>%
    fread(skip = 4,
          fill = TRUE) %>% filter(P < 0.05) %>%
    arrange(P) %>%
    pull(FULL_NAME)
}, simplify = "vector")

intersect_of_significant_gene_sets <-
  Reduce(intersect, gene_set_results)
print(intersect_of_significant_gene_sets)
all_significant_gene_sets <-
  gene_set_results %>% unlist() %>% unique()
print(all_significant_gene_sets)
occurances_of_significant_gene_sets <-
  gene_set_results %>% unlist() %>% table %>% as.data.frame() %>% arrange(desc(Freq))
print(occurances_of_significant_gene_sets)
