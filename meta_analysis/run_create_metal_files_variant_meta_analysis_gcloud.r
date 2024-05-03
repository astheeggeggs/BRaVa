#!/bin/Rscript
library(data.table)
library(dplyr)

# Names of the phenotypes to run
data_dir <- "~/Repositories/BRaVa_curation/data/meta_analysis/gcloud"
out_meta_results_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/meta_results/variant")
system(paste("mkdir -p", out_meta_results_dir))
source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")
source("~/Repositories/BRaVa_curation/phenotypes/BRaVa_phenotypes_utils.r")

# Assumes that we have the files locally within a file structure as defined in the munging scripts
# First, let's determine the collection of phenotypes that we are testing and the 
# collection of (population, phenotype) pairs for that biobank to include in the meta-analysis

biobanks <- dir(data_dir)
results_dt_list <- list()

for (biobank in biobanks)
{
	biobank_results_files <- dir(paste0(data_dir, "/", biobank, "/cleaned/variant"))
	biobank_results_files_full <- dir(paste0(data_dir, "/", biobank, "/cleaned/variant"), full.names=TRUE)
	results <- lapply(biobank_results_files, extract_file_info)
	results_dt_list[[biobank]] <- data.table(filename = biobank_results_files_full, phenotypeID = sapply(results, `[`, 4), pop = sapply(results, `[`, 7))
}

results_dt <- rbindlist(results_dt_list)
pilot_phenotypes <- extract_BRaVa_pilot_phenotypes()
phenotypeIDs <- intersect(pilot_phenotypes, unique(unlist(results_dt$phenotypeID)))

for (phe in phenotypeIDs)
{
	files_variant <- (results_dt %>% filter(phenotypeID == phe))$filename
	files_variant <- paste(files_variant, collapse=",")

	cat(paste0("carrying out meta-analysis of ", phe, "\n"))
	cat(paste0("creating METAL file...\n"))

	system(paste(
		"Rscript create_metal_script.r ",
		"--phenotypeID", phe,
		"--files", files_variant,
		"--out_folder", out_meta_results_dir
		)
	)
	cat(paste0("METAL file creation for meta-analysis of ", phe, " completed\n\n"))
}
