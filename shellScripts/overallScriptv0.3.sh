# Define variables
folderName="als"

# This only needs to be run once
Rscript R/00_install_packages/00_install_packages.R

# This is bespoke to my analysis
Rscript R/convertRawInputsToGeneList/01_convert_raw_inputs_to_gene_list.R

# Gene Set Enrichment using EnrichR
Rscript R/01_gene_set_enrichment/01_gene_enrichment_analysis.R $folderName.csv
Rscript R/01_gene_set_enrichment/02_filter_out_non_significant_go_terms.R $folderName
Rscript R/01_gene_set_enrichment/03_filter_out_redundant_go_terms.R $folderName
Rscript R/01_gene_set_enrichment/04_converting_go_terms_to_gene_lists.R $folderName

# Analysing DGlinker results
Rscript R/02_dglinker_results_analysis/01_analyse_dglinker_gene_lists.R $folderName
Rscript R/02_dglinker_results_analysis/02_extract_best_models_results.R $folderName

# Analysing user input gene lists
Rscript R/03_user_input_gene_lists_analysis/01_convert_gene_sets_to_hcnc_gene_symbols.R $folderName

# Formatting gene sets
Rscript R/04_format_gene_sets/01_identify_disease_related_genes.R $folderName
Rscript R/04_format_gene_sets/02_convert_disease_related_gene_lists_to_entrez_ids.R $folderName
Rscript R/04_format_gene_sets/03_convert_non_disease_related_gene_lists_to_entrez_ids.R $folderName
Rscript R/04_format_gene_sets/04_convert_entrez_id_lists_to_magma_input_including_known_disease_genes.R $folderName
Rscript R/04_format_gene_sets/05_convert_entrez_id_lists_to_magma_input_excluding_known_disease_genes.R $folderName

# Run Magma
mkdir -p data/07_magma_results/$folderName/

software/01_magmaSoftware/magma --annotate window=0 --snp-loc data/01_data_input/07_reference_population/$folderName/g1000_eur.bim  --gene-loc data/01_data_input/08_gene_locations/$folderName/NCBI37.3.gene.loc --out data/07_magma_results/$folderName/gene_annotation

software/01_magmaSoftware/magma --bfile data/01_data_input/07_reference_population/$folderName/g1000_eur --pval data/01_data_input/06_gwas_summary_statistics/$folderName/GCST90027164_buildGRCh37.tsv ncol=N_effective use='rsid,p_value' \
--gene-annot data/07_magma_results/$folderName/gene_annotation.genes.annot --out data/07_magma_results/$folderName/gene_analysis 

directory="data/06_format_gene_sets/"$folderName"/02_magmaInput"

for file in "$directory"/*; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        filename_noext="${filename%.*}"
        software/01_magmaSoftware/magma --gene-results data/07_magma_results/$folderName/gene_analysis.genes.raw \
        --set-annot $file --out data/07_magma_results/$folderName/$filename_noext
    fi
done




sh shellScripts/02_magma/01_magma.sh

# Analyse Magma results
Rscript R/05_magma_results_analysis/ 