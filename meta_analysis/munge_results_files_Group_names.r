#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")

main <- function(args)
{
	folder_to_check <- args$folder
	files_to_check <- grep("\\.gene\\.", dir(folder_to_check, full.names=TRUE), value=TRUE)
	files_to_check <- setdiff(files_to_check, grep("cleaned", files_to_check, value=TRUE))

	for (f in files_to_check)
	{
		cat(paste0("checking file...", f, "...\n"))
		folder <- gsub("(.*\\/)(.*)", "\\1", f)
		file <- gsub("(.*\\/)(.*)", "\\2", f)
		file_info <- extract_file_info(file)
		cat(paste0("file information:\n", paste0(names(file_info), ": ", file_info, collapse="\n"), "\n\n"))
		renamed <- FALSE
		to_check <- c("dataset", "phenotype", "ancestry", "sex")

		for (check in to_check) {
			if (!(file_info[[check]] %in% names(file_check_information[[check]]))) {
				rename_attempt <- TRUE
	            cat(paste0("attempting to rename ", file_info[[check]], "\n"))
	            for (n in names(file_check_information[[check]])) {
					if (file_info[[check]] %in% file_check_information[[check]][[n]]) {
						cat(paste("Renamed:", file_info[[check]], "->", n, "\n"))
						file_info[[check]] <- n
						renamed <- TRUE
						break
					}
				}
			} else {
				rename_attempt <- FALSE
				cat(paste0("No renaming required for ", file_info[[check]], "\n"))
			}

			if (rename_attempt & !renamed) {
				cat(paste0("Cannot find the appropriate renaming for ", file_info[[check]], "\n"))
			}
		}

        if (file_info$type != "gene") {
            cat("This checking script is expecting gene level results\n")
            stop(paste0(file_info$type, " is not 'gene'"))
        }

        cat("\nnaming convention checking complete!\n")

        if (!file_info$gz) {
            cat("\nThe file does not appear to be compressed\n",
                "consider compressing before uploading to the cloud\n")
        }

        # Copy the file to a new location with the new name and run all the tests on the
        # proposed filename.
        if (file_info$binary) {
        	required_info <- c("dataset", "last_name", "analysis_name",
        		"phenotype", "freeze_number",  "sex", "ancestry",
        		"n_cases", "n_controls", "software", "type", "date", "split")
        } else {
        	required_info <- c("dataset", "last_name", "analysis_name",
        		"phenotype", "freeze_number",  "sex", "ancestry",
        		"n", "software", "type", "date", "split")
        }

        # Print the old filename
        if (file_info$gz) {
        	file_info$gz <- "gz"
        	required_info <- c(required_info, "gz")
        }

        new_f <- paste0(folder,
        	paste(file_info[required_info], collapse="."))
        cat(paste0(f, " ->\n ", new_f))

        # Now run all the tests on the proposed file names
		exit_status <- system(paste0("Rscript ~/Repositories/BRaVa_curation/gene_results_file_checks.r --just_filename --file_paths ", new_f))
		if (exit_status == 0) {
			cat("successful, rename\n")
			file.rename(f, new_f)
			cat("continue with further checks...\n")
			dt <- fread(new_f)
			for (n in names(renaming_header_list)) {
				if (!(n %in% names(dt))) {
					cat(paste("attempting to rename to", n,"\n"))
					if (sum(renaming_header_list[[n]] %in% names(dt)) == 1) {
						to_rename <- which(names(dt) %in% renaming_header_list[[n]])
						cat(paste("renamed:", names(dt)[to_rename], "->", n, "\n"))
						names(dt)[to_rename] <- n
					} else {
						cat(paste("cannot find the column to rename to", n, "\n"))
					}
				}
			}

			for (n in names(renaming_group_list)) {
				to_rename <- which(dt$Group %in% renaming_group_list[[n]])
				if (length(to_rename) > 0) {
					cat("Groups (annotation names) to be renamed...\n")
					cat(paste0(dt$Group[to_rename][1], " -> ", n, "\n"))
					dt$Group[to_rename] <- n
				}
			}

			# Check to ensure that all of the gene names are ensembl IDs
			if (all(grepl("ENSG", dt$Region))) {
				cat("all region names contain ensembl gene IDs\n")
				if (all(grepl("_", dt$Region))) {
					cat("looks like the user has combined ensembl gene ID and genesymbol\n")
					dt[, Region:=gsub(".*_(ENSG[0-9]*)$", "\\1", Region)]
					dt[, Region:=gsub("^(ENSG[0-9]*)_.*", "\\1", Region)]
				}
			} else {
				cat("attempting to resolve to ensembl gene ID\n")
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
			type <- ifelse(file_info$binary, "binary", "continuous")

			if (any(!(default_gene_result_columns[[type]] %in% names(dt)))) {
				cat("expected columns are missing!\n")
				print(setdiff(default_gene_result_columns, names(dt)))
				for (colname in setdiff(default_gene_result_columns[[type]], names(dt))) {
					dt[[colname]] <- NA
				}
			} else {
				cat("all required columns are present!\n")
			}

			if (all((unique(dt$Group) %in% correct_names) | is.na(unique(dt$Group)))) {
				cat("all annotations have the correct names...\n")
				out_f <- gsub(folder, paste0(args$out_folder, "/"), new_f)
				fwrite(dt, file=out_f, sep="\t", quote=FALSE)
				cat("file written to:\n")
				cat(paste0(out_f, "\n"))
			}

		} else {
			stop("Error: file name is not compatible with the requested format")
		}
	}
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--folder", default=NULL, required=TRUE,
    help="folder to check")
parser$add_argument("--out_folder", default=NULL, required=TRUE,
    help="folder to export results to")
args <- parser$parse_args()

main(args)

