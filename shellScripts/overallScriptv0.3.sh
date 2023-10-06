# Define variables
folderName="als"

# This only needs to be run once
# Rscript R/00_install_packages/00_install_packages.R

# This is bespoke to my analysis
# Rscript R/convertRawInputsToGeneList/01_convert_raw_inputs_to_gene_list.R

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

# Quality Control GWAS data
Rscript R/05_gwas_data_quality_control/01_gwas_quality_control.R $folderName

# Quality Control Refeence Population Data using plink
files=$(find data/01_data_input/07_reference_population/$folderName/ -maxdepth 1 -type f)
file_names=$(for file in $files; do basename "$file" | sed 's/\..*//'; done)
most_common_name=$(echo "$file_names" | tr ' ' '\n' | sort | uniq -c | sort -nr | head -n 1 | awk '{print $2}')

mkdir -p data/08_quality_controlled_reference_data/$folderName/

/scratch/users/k20064105/statisticalGenetics/software/01_plink/plink \
    --bfile data/01_data_input/07_reference_population/$folderName/$most_common_name \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out data/08_quality_controlled_reference_data/$folderName/ref.QC

/scratch/users/k20064105/statisticalGenetics/software/01_plink/plink \
    --bfile data/01_data_input/07_reference_population/$folderName/$most_common_name \
    --keep data/08_quality_controlled_reference_data/$folderName/ref.QC.fam \
    --extract data/08_quality_controlled_reference_data/$folderName/ref.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out data/08_quality_controlled_reference_data/$folderName/ref.QC

/scratch/users/k20064105/statisticalGenetics/software/01_plink/plink \
    --bfile data/01_data_input/07_reference_population/$folderName/$most_common_name \
    --extract data/08_quality_controlled_reference_data/$folderName/ref.QC.prune.in \
    --keep data/08_quality_controlled_reference_data/$folderName/ref.QC.fam \
    --het \
    --out data/08_quality_controlled_reference_data/$folderName/ref.QC

Rscript R/06_reference_population_quality_control/01_remove_highly_heterozygosity_samples.R $folderName
Rscript R/06_reference_population_quality_control/02_resolve_mismatching_snps.R $folderName

/scratch/users/k20064105/statisticalGenetics/software/01_plink/plink \
    --bfile data/01_data_input/07_reference_population/$folderName/$most_common_name \
    --extract data/08_quality_controlled_reference_data/$folderName/ref.QC.prune.in \
    --keep data/08_quality_controlled_reference_data/$folderName/ref.valid.sample \
    --check-sex \
    --out data/08_quality_controlled_reference_data/$folderName/ref.QC

/scratch/users/k20064105/statisticalGenetics/software/01_plink/plink \
    --bfile data/01_data_input/07_reference_population/$folderName/$most_common_name \
    --extract data/08_quality_controlled_reference_data/$folderName/ref.QC.prune.in \
    --keep data/08_quality_controlled_reference_data/$folderName/ref.QC.valid \
    --rel-cutoff 0.125 \
    --out data/08_quality_controlled_reference_data/$folderName/ref.QC

/scratch/users/k20064105/statisticalGenetics/software/01_plink/plink \
    --bfile data/01_data_input/07_reference_population/$folderName/$most_common_name \
    --make-bed \
    --keep data/08_quality_controlled_reference_data/$folderName/ref.QC.rel.id \
    --out data/08_quality_controlled_reference_data/$folderName/ref.QC \
    --extract data/08_quality_controlled_reference_data/$folderName/ref.QC.snplist \
    --exclude data/08_quality_controlled_reference_data/$folderName/ref.mismatch \
    --a1-allele data/08_quality_controlled_reference_data/$folderName/ref.a1


# Run Magma
mkdir -p data/09_magma_results/$folderName/

software/02_magma/magma --annotate window=0 --snp-loc data/08_quality_controlled_reference_data/$folderName/ref.QC.bim \
--gene-loc $(ls data/01_data_input/08_gene_locations/$folderName/*.gene.loc) --out data/09_magma_results/$folderName/gene_annotation

software/02_magma/magma --bfile data/08_quality_controlled_reference_data/$folderName/ref.QC \
 --pval data/07_quality_controlled_gwas_data/$folderName/quality_controlled_gwas_data.csv ncol=N_effective use='rsid,p_value' \
--gene-annot data/09_magma_results/$folderName/gene_annotation.genes.annot --out data/09_magma_results/$folderName/gene_analysis 

directory="data/06_format_gene_sets/"$folderName"/02_magmaInput"

for file in "$directory"/*; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        filename_noext="${filename%.*}"
        software/02_magma/magma --gene-results data/09_magma_results/$folderName/gene_analysis.genes.raw \
        --set-annot $file --out data/09_magma_results/$folderName/$filename_noext
    fi
done


# Analyse Magma results
#Rscript R/05_magma_results_analysis/ 