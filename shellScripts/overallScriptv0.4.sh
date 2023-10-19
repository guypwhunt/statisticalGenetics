#!/bin/bash
## Define variables
folderName="als"
mafValues=(0.01 0.001 0.0001)

#0.005 0.0025 0.001 0.0005 0.00025 0.0001)

threads=30
numberOfMafValues=${#mafValues[@]}
maxMafProcesses=$((threads<numberOfMafValues?threads:numberOfMafValues))
maxLdscProcesses=$((threads/numberOfMafValues))

## This only needs to be run once
## Rscript R/00_install_packages/00_install_packages.R

## This is bespoke to my analysis
## Rscript R/convertRawInputsToGeneList/01_convert_raw_inputs_to_gene_list.R

## Gene Set Enrichment using EnrichR
Rscript R/01_gene_set_enrichment/01_gene_enrichment_analysis.R $folderName.csv
Rscript R/01_gene_set_enrichment/02_filter_out_non_significant_go_terms.R $folderName
Rscript R/01_gene_set_enrichment/03_filter_out_redundant_go_terms.R $folderName
Rscript R/01_gene_set_enrichment/04_converting_go_terms_to_gene_lists.R $folderName

## Analysing DGlinker results
Rscript R/02_dglinker_results_analysis/01_analyse_dglinker_gene_lists.R $folderName
Rscript R/02_dglinker_results_analysis/02_extract_best_models_results.R $folderName

## Analysing user input gene lists
Rscript R/03_user_input_gene_lists_analysis/01_convert_gene_sets_to_hcnc_gene_symbols.R $folderName

## Formatting gene sets
Rscript R/04_format_gene_sets/01_identify_disease_related_genes.R $folderName
Rscript R/04_format_gene_sets/02_convert_disease_related_gene_lists_to_entrez_ids.R $folderName
Rscript R/04_format_gene_sets/03_convert_non_disease_related_gene_lists_to_entrez_ids.R $folderName
Rscript R/04_format_gene_sets/04_convert_entrez_id_lists_to_magma_input_including_known_disease_genes.R $folderName
Rscript R/04_format_gene_sets/05_convert_entrez_id_lists_to_magma_input_excluding_known_disease_genes.R $folderName
Rscript R/04_format_gene_sets/06_create_gene_coord_file_for_ldsc.R $folderName
Rscript R/04_format_gene_sets/07_create_gene_set_file_for_ldsc_inlcuding_kown_disease_related_genes.R $folderName
Rscript R/04_format_gene_sets/08_create_gene_set_file_for_ldsc_inlcuding_kown_disease_related_genes.R $folderName

qc_gwas_data () {
  Rscript R/05_gwas_data_quality_control/01_gwas_quality_control.R $folderName $maf
}

qc_ref_data () {
files=$(find data/01_data_input/07_reference_population/$folderName/ -maxdepth 1 -type f)
  file_names=$(for file in $files; do basename "$file" | sed 's/\..*//'; done)
  most_common_name=$(echo "$file_names" | tr ' ' '\n' | sort | uniq -c | sort -nr | head -n 1 | awk '{print $2}')

  mkdir -p data/08_quality_controlled_reference_data/$folderName/$maf

  software/01_plink/plink \
      --bfile data/01_data_input/07_reference_population/$folderName/$most_common_name \
      --maf $maf \
      --hwe 1e-6 \
      --geno 0.01 \
      --mind 0.01 \
      --write-snplist \
      --make-just-fam \
      --out data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC

  software/01_plink/plink \
      --bfile data/01_data_input/07_reference_population/$folderName/$most_common_name \
      --keep data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC.fam \
      --extract data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC.snplist \
      --indep-pairwise 200 50 0.25 \
      --out data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC

  software/01_plink/plink \
      --bfile data/01_data_input/07_reference_population/$folderName/$most_common_name \
      --extract data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC.prune.in \
      --keep data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC.fam \
      --het \
      --out data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC

  Rscript R/06_reference_population_quality_control/01_remove_highly_heterozygosity_samples.R $folderName $maf

  software/01_plink/plink \
      --bfile data/01_data_input/07_reference_population/$folderName/$most_common_name \
      --extract data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC.prune.in \
      --keep data/08_quality_controlled_reference_data/$folderName/$maf/ref.valid.sample \
      --check-sex \
      --out data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC

  Rscript R/06_reference_population_quality_control/02_resolve_mismatching_snps.R $folderName $maf

  software/01_plink/plink \
      --bfile data/01_data_input/07_reference_population/$folderName/$most_common_name \
      --extract data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC.prune.in \
      --keep data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC.valid \
      --rel-cutoff 0.125 \
      --out data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC

  Rscript R/06_reference_population_quality_control/03_sample_sex_check.R $folderName $maf

  software/01_plink/plink \
      --bfile data/01_data_input/07_reference_population/$folderName/$most_common_name \
      --make-bed \
      --keep data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC.rel.id \
      --out data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC \
      --extract data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC.snplist \
      --exclude data/08_quality_controlled_reference_data/$folderName/$maf/ref.mismatch \
      --a1-allele data/08_quality_controlled_reference_data/$folderName/$maf/ref.a1
}

magma_analysis () {
mkdir -p data/09_magma_results/$folderName/$maf/

  software/02_magma/magma --annotate window=1,0 --snp-loc data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC.bim \
  --gene-loc $(ls data/01_data_input/08_gene_locations/$folderName/*.gene.loc) --out data/09_magma_results/$folderName/$maf/gene_annotation

  software/02_magma/magma --bfile data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC \
   --pval data/07_quality_controlled_gwas_data/$folderName/$maf/quality_controlled_gwas_data.tsv ncol=N_effective use='rsid,p_value' \
  --gene-annot data/09_magma_results/$folderName/$maf/gene_annotation.genes.annot --out data/09_magma_results/$folderName/$maf/gene_analysis

  directory="data/06_format_gene_sets/"$folderName"/02_magmaInput"

  for file in "$directory"/*; do
    ((n=n%maxLdscProcesses)); ((n++==0)) && wait
      if [ -f "$file" ]; then
          filename=$(basename "$file")
          filename_noext="${filename%.*}"
          software/02_magma/magma --gene-results data/09_magma_results/$folderName/$maf/gene_analysis.genes.raw \
          --set-annot $file --out data/09_magma_results/$folderName/$maf/$filename_noext
      fi &
  done

  Rscript R/07_magma_results_analysis/01_apply_p_value_correction_to_magma_results.R $folderName $maf
  Rscript R/07_magma_results_analysis/02_visulse_magma_results.R $folderName $maf
}

ldsc_analysis () {
  ## Run ldsc
  mkdir -p data/12_ldsc_results/$folderName/$maf/
  Rscript R/08_convert_gwas_data_to_ldsc_format/01_convert_gwas_data_to_ldsc_format.R $folderName $maf

  ## Convert to munge sumstats format
  python software/04_ldsc/munge_sumstats.py \
  --sumstats data/12_ldsc_results/$folderName/$maf/quality_controlled_gwas_data.tsv \
  --out data/12_ldsc_results/$folderName/$maf/quality_controlled_gwas_data \
  --a1-inc \
  --chunksize 500000 \
  --maf-min $maf \
  --info-min 0.8

  ## Calculate unpartitioned LD scores
  python software/04_ldsc/ldsc.py \
	--bfile data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC \
	--l2 \
	--ld-wind-kb 1 \
	--out data/12_ldsc_results/$folderName/$maf/weights

  subFolderNames=("geneSetsExcludingKnownDiseaseGenes" "geneSetsIncludingKnownDiseaseGenes")

  for subFolderName in ${subFolderNames[@]}
  do
    outputDirectory=data/12_ldsc_results/$folderName/$maf/$subFolderName/
    parentDirectory=data/06_format_gene_sets/$folderName/03_ldscInput/$subFolderName/

    directories=$(ls $parentDirectory)
    for directory in ${directories[@]}
    do
      subOutputDirectory=$outputDirectory$directory
      mkdir -p $subOutputDirectory

      subdirectory=$parentDirectory$directory/
      fileNames=$(ls $subdirectory)
      for fileName in ${fileNames[@]}
      do
      ((x=x%maxLdscProcesses)); ((x++==0)) && wait
      filePath=$subdirectory$fileName

      outputFileName="${fileName%.*}"

      python software/04_ldsc/make_annot.py \
  		--gene-set-file $filePath \
  		--gene-coord-file data/06_format_gene_sets/$folderName/03_ldscInput/geneCoordFile.txt \
  		--windowsize 1000 \
  		--bimfile data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC.bim \
  		--annot-file $subOutputDirectory/$outputFileName.annot.gz

  		python software/04_ldsc/ldsc.py \
  		--l2 \
  		--bfile data/08_quality_controlled_reference_data/$folderName/$maf/ref.QC \
  		--ld-wind-kb 1 \
  		--annot $subOutputDirectory/$outputFileName.annot.gz \
  		--thin-annot \
  		--out $subOutputDirectory/$outputFileName

  		python software/04_ldsc/ldsc.py \
    	--h2 data/12_ldsc_results/$folderName/$maf/quality_controlled_gwas_data.sumstats.gz \
    	--ref-ld $subOutputDirectory/$outputFileName \
    	--w-ld data/12_ldsc_results/$folderName/$maf/weights \
    	--overlap-annot \
    	--not-M-5-50 \
    	--print-coefficients \
    	--out $subOutputDirectory/$outputFileName &
      done
    done
  done

  Rscript R/09_ldsc_results_analysis/01_consolidate_ldsc_results.R $folderName $maf
  Rscript R/09_ldsc_results_analysis/02_apply_p_value_correction_to_ldsc_results.R $folderName $maf
  Rscript R/09_ldsc_results_analysis/03_visulse_ldsc_results.R $folderName $maf


}

prset_analysis () {
  mkdir -p data/11_prsice_results/$folderName/$maf
}

process_maf () {
  qc_gwas_data $maf
  qc_ref_data $maf
  magma_analysis $maf
  ldsc_analysis $maf
  prset_analysis $maf
}

for maf in "${mafValues[@]}"; do
  ((i=i%maxMafProcesses)); ((i++==0)) && wait
  process_maf $maf &
done

