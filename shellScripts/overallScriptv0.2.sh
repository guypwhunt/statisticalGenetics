# This only needs to be run once
Rscript R/00_install_packages/00_install_packages.R

# This is bespoke to my analysis
Rscript R/convertRawInputsToGeneList/01_convert_raw_inputs_to_gene_list.R

# Gene Set Enrichment
Rscript R/01_gene_set_enrichment/01_gene_enrichment_analysis.R als.csv
Rscript R/01_gene_set_enrichment/02_filter_out_non_significant_go_terms.R als
Rscript R/01_gene_set_enrichment/03_filter_out_redundant_go_terms.R als
Rscript R/01_gene_set_enrichment/04_converting_go_terms_to_gene_lists.R als


sh shellScripts/02_magma/01_magma.sh
