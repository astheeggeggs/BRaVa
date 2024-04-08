#!/bin/Rscript

# Names of the phenotypes to run
phenotypes <- c("Coronary_artery_disease")
biobanks <- c("uk-biobank", "genes-and-health", "mgbb", "gel", "pmbb", "bbj", "all-of-us")
system("mkdir -p ../data/meta_analysis/gcloud")
data_dir <- "../data/meta_analysis/gcloud"
out_meta_results_dir <- "../data/meta_analysis/gcloud/meta_results"
system(paste("mkdir -p", out_meta_results_dir))
download <- FALSE

for (phenotype in phenotypes)
{
	if (download) {
		for (biobank in biobanks) {
			# Download the required data for the meta-analysis
			system(paste0(
				"gsutil cp gs://brava-meta-upload-", biobank, "/*", phenotype, "* ",
				data_dir, "/")
			)
			system(paste0(
				"gsutil cp gs://brava-meta-upload-", biobank, "/pilot/*", phenotype, "* ",
				data_dir, "/")
			)
		}
	}

	system("Rscript munge_results_files_Group_names.r")

	files <- grep(phenotype, dir(data_dir, full.names=TRUE), value=TRUE)
	files_gene <- paste(grep("gene\\..*cleaned.*", files, value=TRUE), collapse=",")
	print(strsplit(files_gene, split=",")[[1]])
	system(paste(
		"Rscript meta_analysis.r --file_paths",
		files_gene, "--out",
		paste0(out_meta_results_dir, "/n_cases_100/", phenotype, "_gene_meta_analysis_100_cutoff.tsv.gz"))
	)
	system(paste(
		"Rscript meta_analysis.r --file_paths",
		files_gene, "--case_control_threshold 1000", "--out",
		paste0(out_meta_results_dir, "/n_cases_1000/", phenotype, "_gene_meta_analysis_1000_cutoff.tsv.gz"))
	)
}
