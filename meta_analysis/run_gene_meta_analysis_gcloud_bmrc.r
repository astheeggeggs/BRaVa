#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

main <- function(args)
{
	data_dir <- args$data_dir
	n_cases <- args$n_cases
	out_dir <- args$out_dir
	phe <- args$phenotypeID

	out_meta_results_dir <- paste0(out_dir, "/gene/n_cases_", n_cases)
	
	# Ensure that the folder is already present
	system(paste("mkdir -p", out_meta_results_dir))
	source("meta_analysis_utils.r")
	source("../phenotypes/BRaVa_phenotypes_utils.r")

	# Assumes that we have the files locally within a file structure as defined in the munging scripts
	# First, let's determine the collection of phenotypes that we are testing and the 
	# collection of (population, phenotype) pairs for that biobank to include in the meta-analysis

	biobanks <- dir(data_dir)
	results_dt_list <- list()

	for (biobank in biobanks)
	{
		biobank_results_files <- dir(paste0(data_dir, "/", biobank, "/cleaned/gene"))
		biobank_results_files_full <- dir(paste0(data_dir, "/", biobank, "/cleaned/gene"), full.names=TRUE)
		results <- lapply(biobank_results_files, extract_file_info)
		results_dt_list[[biobank]] <- data.table(filename = biobank_results_files_full, phenotypeID = sapply(results, `[[`, 4), pop = sapply(results, `[[`, 7))
	}

	results_dt <- rbindlist(results_dt_list)

	if (is.null(phe)) {
		phes <- c(
			"AMD", "Asth", "AFib", "BenCervUterNeo", "BenIntNeo",
			"BreastCanc", "CervCanc", "COPD", "CRF", "ColonRectCanc",
			"CAD", "EFRMB", "FemInf", "Gout", "HF", "HTN", "IBD",
			"IFHern", "ILDSarc", "MatHem", "NonRheuValv", "Pancreat",
			"PeptUlcer", "PAD", "Psori", "RheumHeaDis", "RheumArth",
			"Stroke", "T2Diab", "Urolith", "VaricVeins", "VTE", "ALT",
			"AlcCons", "AST", "BMI", "CRP", "HDLC", "Height", "LDLC",
			"TChol", "TG", "WHRBMI", "HipRep"
		)
	} else {
		phes <- phe
	}

	for (phe in phes) {
		files_gene <- (results_dt %>% filter(phenotypeID == phe))$filename
		files_gene <- paste(files_gene, collapse=",")
		out <- paste0(out_meta_results_dir, "/", phe, "_gene_meta_analysis_", n_cases, "_cutoff.tsv.gz")
		cat(paste0("carrying out meta-analysis of ", phe, "\n"))
		cat(paste0("\nFiles in the analysis: ",
			paste0(strsplit(files_gene, split=",")[[1]], collapse='\n'), "\n"))
		system(paste(
			"sbatch run_meta_analysis_gcloud_bmrc.sh",
			files_gene, out))
		cat(paste0("submitted meta-analysis of ", phe, " completed\n\n"))
	}
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--data_dir", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs",
	required=FALSE, help="Location of the folder containing the input into the meta-analysis")
parser$add_argument("--n_cases", default=100, required=FALSE,
	help="Minimum number of cases")
parser$add_argument("--out_dir", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs",
	required=FALSE, help="Output folder path")
parser$add_argument("--phenotypeID", required=FALSE, default=NULL,
	help="The phenotype ID to run meta-analysis on. Note: thus exactly must match the naming in input folder.")
args <- parser$parse_args()

main(args)