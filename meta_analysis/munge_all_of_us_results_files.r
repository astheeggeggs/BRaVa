#!/bin/Rscript
library(data.table)
library(dplyr)

# For naming files
biobank <- "all-of-us"
last_name <- "Lu"
analysis_name <- "pilot"
freeze_number <- "JULY23Freeze"

data_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/", biobank, "/raw")
out_data_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/", biobank, "/cleaned")
system(paste0("mkdir -p ", data_dir, "/gene"))
system(paste0("mkdir -p ", data_dir, "/variant"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- TRUE

data_dictionary <- fread("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/all-of-us/brava_44_pilot_phenotype_filename.csv")

files_to_retain <- data.table(data_dictionary %>% 
	group_by(brava_code, pheno_sex, pop) %>% 
	filter(n_sample_label == max(n_sample_label)))

# Note that the following files are the same (as ascertained by md5sum):
# The reason for this is likely that Wenhan uploaded some files and then after the
# fact, noticed that our requested phenotype labels didn't quite match.

# Total cholesterol
# all-of-us.Lu.pilot.cholesterol.ALL.{AFR,AMR,EAS,EUR,SAS}.{count}.SAIGE.gene.20240401.txt.gz
# all-of-us.Lu.pilot.TChol_3027114.ALL.{AFR,AMR,EAS,EUR,SAS}.{count}.SAIGE.gene.20240401.txt.gz

# HDL cholesterol
# all-of-us.Lu.pilot.HDL.ALL.{AFR,AMR,EAS,EUR,SAS}.{count}.SAIGE.gene.20240401.txt.gz
# all-of-us.Lu.pilot.HDLC_3007070.ALL.{AFR,AMR,EAS,EUR,SAS}.{count}.SAIGE.gene.20240401.txt.gz

# LDL cholesterol
# all-of-us.Lu.pilot.LDL.ALL.{AFR,AMR,EAS,EUR,SAS}.{count}.SAIGE.gene.20240401.txt.gz
# all-of-us.Lu.pilot.LDLC_3028288.ALL.{AFR,AMR,EAS,EUR,SAS}.{count}.SAIGE.gene.20240401.txt.gz

# T2D phenotypes
# all-of-us.Lu.pilot.T2Diab_EM_202.2.ALL.{AFR,AMR,EAS,EUR,SAS}.{case_count}.{control_count}.SAIGE.gene.20240424.txt.gz
# all-of-us.Lu.pilot.T2D_EM_202.2.ALL.{AFR,AMR,EAS,EUR,SAS}.{case_count}.{control_count}.SAIGE.gene.20240424.txt.gz
# all-of-us.Lu.pilot.T2Diab_250.2.ALL.{AFR,AMR,EAS,EUR,SAS}.{case_count}.{control_count}.SAIGE.gene.20240424.txt.gz
# all-of-us.Lu.pilot.T2D_250.2.ALL.{AFR,AMR,EAS,EUR,SAS}.{case_count}.{control_count}.SAIGE.gene.20240424.txt.gz

if (download) {
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank, "/*.gene.* ",
		data_dir, "/gene/")
	)
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank, "/*.variant.* ",
		data_dir, "/variant/")
	)
}

files_to_retain <- files_to_retain %>% 
	mutate(new_gene_file_name = paste(
		biobank,
		last_name,
		analysis_name,
		brava_code,
		freeze_number,
		pheno_sex,
		toupper(pop),
		ifelse(trait_type=="continuous", n_cases, paste(n_cases, n_controls, sep=".")),
		"SAIGE",
		"gene",
		gsub("^.*gene.([0-9]+.*$)", "\\1", gene_file_name),
		sep=".")) %>% 
	mutate(new_variant_file_name = gsub(
		"\\.gene\\.", "\\.variant\\.", new_gene_file_name))

from <- paste0(data_dir, "/gene/", files_to_retain$gene_file_name)
to <- paste0(out_data_dir, "/gene/", files_to_retain$new_gene_file_name)

file.rename(from, to)
system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", out_data_dir, "/gene",
	" --type ", "gene",
	" --write",
	" --out_folder ", out_data_dir, "/gene")
)

from <- paste0(data_dir, "/variant/", files_to_retain$var_file_name)
to <- paste0(out_data_dir, "/variant/", files_to_retain$new_variant_file_name)

file.rename(from, to)
system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", out_data_dir, "/variant",
	" --type ", "variant",
	" --write",
	" --out_folder ", out_data_dir, "/variant")
)

