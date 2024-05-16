#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

main <- function(args)
{
	data_dir <- args$analysis_results_folder
	n_cases <- args$n_cases
	out_dir <- args$out_dir
	phe <- args$phenotypeID

	out_plot_dir <- paste0(out_dir, "/gene/n_cases_", n_cases)
	
	# Ensure that the folder is already present
	system(paste("mkdir -p", out_plot_dir))
	source("meta_analysis_utils.r")
	source("../phenotypes/BRaVa_phenotypes_utils.r")

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

	datasets <- c("all-of-us", "alspac", "biome", "bbj", "ckb", "ccpm",
		"decode", "egcut", "dan-rav", "genes-and-health", "gel", "pmbb",
		"mgbb", "qatar-genomes", "uk-biobank", "viking-genes")

	ancestries <- c("AFR", "AMR", "EAS", "EUR", "SAS")

	for (phe in phes) {
		for (dataset in datasets) {
			for (anc in ancestries) {
				file_gene <- grep(
					paste0(".*cleaned.*", dataset, "\\..*",  phe, ".*", anc),
					dir(data_dir, full.names=TRUE, recursive=TRUE), value=TRUE)
				out <- paste0(out_plot_dir, "/", dataset, "_", phe, "_", anc, "_gene_meta_analysis_qq.pdf")
				cat(paste0("carrying out plotting of gene QQ for ", phe, "\n"))
				if (length(file_gene == 1)) {
					cat(paste0("using file: ", file_gene, "\n"))
					system(paste(
						"sbatch run_analysis_qq_gcloud_bmrc.sh",
						file_gene, out))
					cat(paste0("submitted meta-analysis QQ plotting of ", phe, "\n\n"))
				}
				if (length(file_gene) > 1) {
					stop("Multiple matches to this (phenotyes, ancestry, biobank) tuple")
				}
			}
		}
	}	
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--analysis_results_folder", required=FALSE,
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs")
parser$add_argument("--out_dir",
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plots_biobank_specific",
	required=FALSE, help="Output folder path")
parser$add_argument("--n_cases", default=100, required=FALSE,
	help="Minimum number of cases")
parser$add_argument("--phenotypeID", required=FALSE, default=NULL,
	help="The phenotype ID to plot. If null, this script will plot everything in the folder.")
args <- parser$parse_args()

main(args)
