library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)

input_file_name <- commandArgs(trailingOnly = TRUE)[1]
mfa <- commandArgs(trailingOnly = TRUE)[2] %>% as.numeric() %>%
  format(scientific = FALSE)

# mfa <- 0.0025
#
# input_file_name <- "als"

input_file_path <-
  paste0("data/14_p_value_adjusted_ldsc_results/",
         input_file_name,
         "/",
         mfa,
         "/")

output_file_path <-
  paste0("data/14_p_value_adjusted_ldsc_results/",
         input_file_name,
         "/",
         mfa,
         "/")

results_including_gwas_hits <-
  fread(paste0(
    input_file_path,
    "geneSetsIncludingKnownDiseaseGenes.csv"
  ))

results_excluding_gwas_hits <-
  fread(paste0(
    input_file_path,
    "geneSetsExcludingKnownDiseaseGenes.csv"
  )) %>%
  select(c(gene_set, Enrichment_p, Adjust_P)) %>%
  rename(ex_p = Enrichment_p, ex_Adjust_P = Adjust_P)


all_results <-
  merge(results_including_gwas_hits,
        results_excluding_gwas_hits,
        by = "gene_set") %>%
  filter(Enrichment_p < 0.05) %>%
  mutate(
    Log_Adjust_P = 0 - log10(Adjust_P),
    gene_set  = str_replace(str_replace_all(gene_set , "_", " "), ".csv", "")
  ) %>%
  arrange(Adjust_P)

insert_line_break <- function(input_string) {
  words <- strsplit(input_string, " ")[[1]]
  if (length(words) >= 4) {
    words[4] <- paste0(words[4], "\n")
  }
  return(paste(words, collapse = " "))
}

all_results$gene_set <-
  lapply(all_results$gene_set, insert_line_break) %>%
  unlist()

all_results$excluding_significance <-
  ifelse(
    all_results$ex_Adjust_P < 0.05,
    "Adjusted P-value Significant",
    ifelse(
      all_results$ex_P < 0.05,
      "Nominal P-value Significant",
      "Not Significant"
    )
  ) %>%
  factor(
    levels = c(
      "Adjusted P-value Significant",
      "Nominal P-value Significant",
      "Not Significant"
    )
  )

ggplot(all_results,
       aes(
         x = reorder(gene_set,-Log_Adjust_P),
         y = Log_Adjust_P,
         color = excluding_significance
       )) +
  geom_point() +
  labs(x = "Gene Set", y = "0-log10(Adjusted P-Value)", size = "Number of Genes") +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    size = 6
  )) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dotted",
    color = "black"
  ) +
  guides(color = guide_legend(title = "Significance excluding\nGWAS and input genes")) +
  ylim(0, NA) +
  scale_color_manual(values = c("Adjusted P-value Significant" = "red",
                                "Nominal P-value Significant"="green",
                                "Not Significant"="blue"))

ggsave(
  paste0(output_file_path, "ldsc_results.png"),
  dpi = 600,
  width = 20,
  height = 15,
  units = "cm"
)
