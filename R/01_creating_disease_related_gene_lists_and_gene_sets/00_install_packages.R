install.packages(c("data.table",
                   "tidyverse",
                   "dplyr",
                   "stringr",
                   "readxl",
                   "ggvenn"),
                 dependencies = TRUE,
                 repos='http://cran.uk.r-project.org')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("biomaRt",
                       "enrichR",
                       "GO.db",
                       "rrvgo",
                       "org.Hs.eg.db"),
                     ask = FALSE)
