#!/bin/Rscript
library(data.table)
library(dplyr)

# Renaming errors
renaming_group_list <- list(
	`damaging_missense_or_protein_altering` = c("damaging_missense", "missenseLC"),
	`other_missense_or_protein_altering` = c("other_missense"),
	`pLoF;damaging_missense_or_protein_altering` = c("pLoF;damaging_missense", "pLoF;missenseLC", "damaging_missense;pLoF"),
	`pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous` = c(
		"pLoF;damaging_missense;other_missense;synonymous",
		"pLoF;damaging_missense;other_missense;synonymous;pLoF"),
	Cauchy = c(NA)
)

renaming_header_list <- list(
	`Region` = c("gene"),
	`Group` = c("annot"),
	`max_MAF` = c("max_maf"),
	`Pvalue` = c("p_value"),
	`Pvalue_Burden` = c("p_value_burden"),
	`Pvalue_SKAT` = c("p_value_skat"),
	`BETA_Burden` = c("beta_burden"),
	`SE_Burden` = c("se_burden"),
	`MAC` = c("mac"),
	`MAC_case` = c("mac_case"),
	`MAC_control` = c("mac_control"), 
	`Number_rare` = c("rare_var_count"),
	`Number_ultra_rare` = c("ultrarare_var_count")
)

correct_names <- c("pLoF",
	"damaging_missense_or_protein_altering",
	"other_missense_or_protein_altering",
	"synonymous",
	"pLoF;damaging_missense_or_protein_altering",
	"pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
	"Cauchy")

folder_to_check <- "~/Repositories/BRaVa_curation/data/meta_analysis/gcloud"
files_to_check <- grep("\\.gene", dir(folder_to_check, full.names=TRUE), value=TRUE)
files_to_check <- setdiff(files_to_check, grep("cleaned", files_to_check, value=TRUE))
print(files_to_check)

# Expected columns
default_gene_result_columns <- c(
	"Region", "Group", "max_MAF", "Pvalue",
	"Pvalue_Burden", "Pvalue_SKAT","BETA_Burden", "SE_Burden",
	"MAC", "MAC_case", "MAC_control", "Number_rare", "Number_ultra_rare")

for (f in files_to_check) {
	cat(paste0("checking file...", f, "...\n"))
	dt <- fread(f)

	for (n in names(renaming_header_list)) {
		if (!(n %in% names(dt))) {
			cat(paste("Attempting to rename to", n,"\n"))
			if (sum(renaming_header_list[[n]] %in% names(dt)) == 1) {
				to_rename <- which(names(dt) %in% renaming_header_list[[n]])
				cat(paste("Renamed:", names(dt)[to_rename], "->", n, "\n"))
				names(dt)[to_rename] <- n
			} else {
				cat(paste("Cannot find the column to rename to", n, "\n"))
			}
		}
	}

	for (n in names(renaming_group_list)) {
		to_rename <- which(dt$Group %in% renaming_group_list[[n]])
		if (length(to_rename) > 0) {
			dt$Group[to_rename] <- n
		}
	}

	# Check to ensure that all of the gene names are ensembl IDs
	if (all(grepl("ENSG", dt$Region))) {
		cat("All Region names contain ensembl gene IDs\n")
		if (all(grepl("_", dt$Region))) {
			cat("Looks like the user has combined ensembl gene ID and genesymbol\n")
			dt[, Region:=gsub(".*_(ENSG[0-9]*)$", "\\1", Region)]
			dt[, Region:=gsub("^(ENSG[0-9]*)_.*", "\\1", Region)]
		}
	} else {
		cat("Attempting to resolve to ensembl gene ID\n")
		dt_hgnc <- fread("230117_hgncid_ensembl.txt.gz", select = c("ensembl_gene_id", "hgnc_symbol"))
		dt_hgnc[, Region:=hgnc_symbol]
		dt_hgnc$Region <- ifelse(dt_hgnc$Region == "", dt_hgnc$ensembl_gene_id, dt_hgnc$Region)
		dt_hgnc[, hgnc_symbol:=NULL]
		setkey(dt_hgnc, "Region")
		setkey(dt)
		dt <- merge(dt, dt_hgnc, all.x=TRUE)
		dt[, Region:=ensembl_gene_id]
		dt <- dt %>% filter(Region != "")
	}

	# Remove all other columns
	dt <- dt %>% select(intersect(names(dt), names(renaming_header_list)))

	if (any(!(default_gene_result_columns %in% names(dt)))) {
		cat("Expected columns are missing!\n")
		print(setdiff(default_gene_result_columns, names(dt)))
		for (colname in setdiff(default_gene_result_columns, names(dt))) {
			dt[[colname]] <- NA
		}
	}
	if (all((unique(dt$Group) %in% correct_names) | is.na(unique(dt$Group)))) {
		fwrite(dt, file=gsub(".txt.gz", ".cleaned.txt.gz", f), sep="\t", quote=FALSE)
	}
}
