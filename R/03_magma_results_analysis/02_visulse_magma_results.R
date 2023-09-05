# Load required packages
library(dplyr)
library(ggplot2)

# Read CSV as a data frame
data <- read.csv("data/02_magma/09_pValueAdjustedGeneSetAnalysis/allGwasAssociatedGenesRemoved/0.csv")

# Filter data based on P values and create Log_Adjust_P column
filtered_data <- data %>%
  filter(P < 0.1) %>%
  mutate(Log_Adjust_P = 0 - log10(Adjust_P))

ggplot(filtered_data, aes(x = reorder(FULL_NAME, -Log_Adjust_P), y = Log_Adjust_P, size = NGENES,  color = BETA)) +
  geom_point() +
  labs(x = "Gene Set", y = "0-log10(Adjusted P-Value)", size = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
  scale_colour_continuous(type = "viridis") +
  coord_cartesian(ylim = c(0, 1.75))

