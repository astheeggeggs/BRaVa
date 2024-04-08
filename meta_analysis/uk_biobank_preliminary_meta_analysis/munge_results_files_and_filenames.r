#!/bin/Rscript

# Code to combine results files across chromosomes and rename files according to the pilot analysis plan

# First, grab the files of interest
phenotypes <- "Coronary_artery_disease"
output_folder <- "~/Repositories/BRaVa_curation/data/sept2023"
pops <- c("AFR", "AMR", "EAS", "EUR", "SAS")
# Determine the date
date <- format(Sys.time(), "%Y%m%d")
date <- "20240110"
dataset <- "uk-biobank"
last_name <- "palmer"
analysis_name <- "PRELIMINARY"
freeze_number <- "JULY23Freeze"
sex <- "ALL"
method <- "SAIGE"

gcloud_bucket <- paste0("brava-meta-upload-", dataset)

phenotype_mapping <- list(
	Coronary_artery_disease = "Coronary_artery_disease")

determine_binary_filename <- function(dataset, last_name, analysis_name, phenotype, sex,
	ancestry, n_cases, n_controls, type, date,
	method="SAIGE", freeze_number="JULY23Freeze") {
	# Format
	# [dataset].[last_name].[analysis_name].[phenotype].[freeze_number].[sex].[ancestry].[n_cases].[n_controls].[SAIGE].{gene,variant}.[YYYYMMDD].txt.gz
	return(paste(
		dataset,
		last_name,
		analysis_name,
		phenotype,
		freeze_number,
		sex,
		ancestry,
		n_cases,
		n_controls,
		method,
		type,
		date,
		"txt.gz",
		sep="."))
}

determine_cts_filename <- function(dataset, last_name, analysis_name, phenotype, sex,
	ancestry, n, type, date,
	method="SAIGE", freeze_number="JULY23Freeze") {
	# Format
	# [dataset].[last_name].[analysis_name].[phenotype].[freeze_number].[sex].[ancestry].[n].[SAIGE].{gene,variant}.[YYYYMMDD].txt.gz
	return(paste(
		dataset,
		last_name,
		analysis_name,
		phenotype,
		freeze_number,
		sex,
		ancestry,
		n,
		method,
		type,
		date,
		"txt.gz",
		sep="."))
}

system(paste("mkdir", output_folder))
for (phenotype in phenotypes) {
	outputs <- system("dx ls brava/outputs/step2/sept2023", intern=TRUE)
	grep(phenotype, outputs, value=TRUE)
	system(paste0("dx download ", "brava/outputs/step2/sept2023/*", phenotype, "* -o ", output_folder))
}

# Make gene and variant level folders
gene_folder <- paste0(output_folder, "/gene")
variant_folder <- paste0(output_folder, "/variant") 
system(paste("mkdir", gene_folder))
system(paste("mkdir", variant_folder))

# Move the variant level files to the variant folder
system(paste0("mv ", output_folder, "/*singleAssoc* ", variant_folder))
# Move the remainder to the gene level folder
system(paste0("mv ", output_folder, "/*gz ", gene_folder))

# Create file for the final output
system(paste("mkdir", paste0(gene_folder, "/combined")))
system(paste("mkdir", paste0(variant_folder, "/combined")))

variant_files <- dir(variant_folder, full.names=TRUE)
gene_files <- dir(gene_folder, full.names=TRUE)

# Combine the variant based files across the chromosomes, for each population label
for (phenotype in phenotypes) {
	for (pop in pops) {
		pop_variant_files <- grep(paste0("_", pop), variant_files, value=TRUE)
		if (length(pop_variant_files) > 0)
		{
			dt_list <- list()
			for (file in pop_variant_files) {
				dt_list[[file]] <- fread(file)
			}
			dt_variant <- rbindlist(dt_list)

			# Determine whether the trait is continuous or binary
			binary <- ifelse("N_case" %in% names(dt_variant), TRUE, FALSE)
			# Determine whether the sex is ALL, M, or F
			sex <- ifelse(grepl("_F.txt", pop_variant_files[1]), "F",
					ifelse(grepl("_M.txt", pop_variant_files[1]), "M", "ALL"))

			pop_gene_files <- grep(paste0("_", pop), gene_files, value=TRUE)
			dt_list <- list()
			for (file in pop_gene_files) {
				dt_list[[file]] <- fread(file)
			}
			dt_gene <- rbindlist(dt_list)

			# Extract the relevant information from the resultant files and rename the files 
			# according to the pilot analysis plan
			if (binary) {
				# Extract the case and control count
				n_cases <- mean(dt_variant$N_case)
				n_controls <- mean(dt_variant$N_ctrl)
				filename_variant <- determine_binary_filename(
					dataset, last_name, analysis_name, phenotype_mapping[[phenotype]],
					sex, pop, n_cases, n_controls, "variant", date)
				filename_gene <- gsub("variant", "gene", filename_variant)
				filename_variant <- paste0(paste0(variant_folder, "/combined/"), filename_variant)
				filename_gene <- paste0(paste0(gene_folder, "/combined/"), filename_gene)
				fwrite(dt_variant, quote=FALSE, file=filename_variant, sep="\t")
				fwrite(dt_gene, quote=FALSE, file=filename_gene, sep="\t")
			} else {
				# Extract n
				n <- mean(dt_variant$N)
				filename_variant <- determine_cts_filename(
					dataset, last_name, analysis_name, phenotype_mapping[[phenotype]],
					sex, pop, n, "variant", date)
				filename_gene <- gsub("variant", "gene", filename_variant)
				filename_variant <- paste0(paste0(variant_folder, "/combined/"), filename_variant)
				filename_gene <- paste0(paste0(gene_folder, "/combined/"), filename_gene)
				fwrite(dt_variant, quote=FALSE, file=filename_variant, sep="\t")
				fwrite(dt_gene, quote=FALSE, file=filename_gene, sep="\t")
			}
		}
	}
}

# Upload the results to the relevant bucket
system(paste0("gsutil cp ", variant_folder, "/combined/* gs://", gcloud_bucket, "/"))
system(paste0("gsutil cp ", gene_folder, "/combined/* gs://", gcloud_bucket, "/"))


