#!/bin/Rscript

# Names of the phenotypes to run
data_dir <- "~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/genes-and-health/raw"
out_data_dir <- "~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/genes-and-health/cleaned"
system(paste0("mkdir -p ", data_dir, "/gene"))
system(paste0("mkdir -p ", data_dir, "/variant"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- FALSE
biobank <- "genes-and-health"
freeze_to_match <- "JULY23Freeze"
round_to_match <- "Round1"

rename_files <- function(results_type="gene")
{
	# Change the filenames with gsub
	from <- list.files(paste0(data_dir, "/", results_type), pattern=freeze_to_match, full.names=TRUE)
	from <- from[which(grepl(paste0(round_to_match, ".*", freeze_to_match), from))]
	phenotype <- gsub(paste0(".*", round_to_match, "\\.(.*)\\.", freeze_to_match, ".*"), "\\1", from)
	new_phenotype <- gsub("\\.", "_", phenotype)
	new_phenotype <- gsub("[_]+", "_", new_phenotype)
	new_phenotype <- gsub("_$", "", new_phenotype)
	start <- gsub(paste0("(.*", round_to_match, "\\.).*"), "\\1", from)
	end <- gsub(paste0(".*(\\.", freeze_to_match, ".*)"), "\\1", from)

	file.rename(from, paste0(start, new_phenotype, end))
}

if (download) {
	system(paste0(
		"gsutil cp gs://brava-meta-upload-", biobank, "/*.gene.* ",
		data_dir, "/gene/")
	)
	rename_files()
	
	system(paste0(
		"gsutil cp gs://brava-meta-upload-", biobank, "/*.variant.* ",
		data_dir, "/variant/")
	)
	rename_files(results_type="variant")
}

system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", data_dir, "/gene",
	" --type ", "gene",
	" --write",
	" --out_folder ", out_data_dir, "/gene")
)

system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", data_dir, "/variant",
	" --type ", "variant",
	" --write",
	" --out_folder ", out_data_dir, "/variant")
)
