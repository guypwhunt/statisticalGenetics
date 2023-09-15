There is no INFO OR MAF column (but there is effect allele frequency column which could be used instead of MAD

There are no duplicated rsIDs



library(data.table)
# Read in file
dat <- fread("/scratch/users/k20064105/statisticalGenetics/data/02_magma/04_gwasSummaryStatistics/GCST90027164_buildGRCh37.tsv")
# Filter out SNPs
result <- dat[effect_allele_frequency > 0.01]
# Output the gz file
fwrite(result, "Height.gz", sep="\t")