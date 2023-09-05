data/02_magma/01_magmaSoftware/magma --annotate window=0 --snp-loc data/02_magma/02_referencePopulation/g1000_eur/g1000_eur.bim  --gene-loc data/02_magma/03_geneLocations/NCBI37.3/NCBI37.3.gene.loc --out data/02_magma/06_annotation/0

data/02_magma/01_magmaSoftware/magma --bfile data/02_magma/02_referencePopulation/g1000_eur/g1000_eur --pval data/02_magma/04_gwasSummaryStatistics/GCST90027164_buildGRCh37.tsv ncol=N_effective use='rsid,p_value' --gene-annot data/02_magma/06_annotation/0.genes.annot --out data/02_magma/07_geneAnalysis/0

data/02_magma/01_magmaSoftware/magma --gene-results data/02_magma/07_geneAnalysis/0.genes.raw --set-annot data/01_geneLists/09_magmaInput/magmaInput.txt --out data/02_magma/08_geneSetAnalysis/0