#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(argparse)

# Pass a comma separated list of filenames
tests <- c("Burden", "SKAT", "SKAT-O")
source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")

# Extract the files to read in
# Pass the name of the trait to be used in the output
# Let the case control threshold be an option
# Let the N threshold be an option

checks <- function(file_info, file_info_template, datasets) {
    if (file_info$phenotype != file_info_template$phenotype) {
        print(file_info)
        stop("phenotype does not match - check files or rename")
    }
    if (file_info$sex != file_info_template$sex) {
        print(file_info)
        stop("meta-analysis of different sex - check files")
    }
    if (file_info$type != file_info_template$type) {
        print(file_info)
        stop("attempting to meta-analyse gene based results with variant results - check files")
    }
    if (file_info$binary != file_info_template$binary) {
        print(file_info)
        stop("attempting to analyse binary phenotype with a continuous phenotype")
    }
}

main <- function(args)
{
    files <- strsplit(args$file_paths, split=",")[[1]]
    case_control_threshold <- as.integer(args$case_control_threshold)
    dt_list <- list()
    folder <- gsub("(.*\\/)(.*)", "\\1", files[1])
    file <- gsub("(.*\\/)(.*)", "\\2", files[1])
    file_info_template <- extract_file_info(file)
    print(args)

    for (f in files)
    {
        print(f)
        # Throw an error if the file is not there
        if (!file.exists(f)) {
            stop(paste("File:", file, "does not exist."))
        }

        # Get information about the file
        folder <- gsub("(.*\\/)(.*)", "\\1", f)
        file <- gsub("(.*\\/)(.*)", "\\2", f)
        file_info <- extract_file_info(file)
        checks(file_info, file_info_template)

        # Throw an error if the case count is too low
        if (file_info$binary) {
            file_info$n_cases <- as.integer(file_info$n_cases)
            file_info$n_controls <- as.integer(file_info$n_controls)
            if (file_info$n_cases < case_control_threshold) {
                next
            }
            if (file_info$n_controls < case_control_threshold) {
                next
            }
        }

        dt_list[[file]] <- fread(cmd = ifelse(file_info$gz,
            paste0("gzcat ", folder, file),
            paste0("cat ", folder, file)))
        dt_list[[file]]$dataset <- file_info$dataset
        dt_list[[file]]$ancestry <- file_info$ancestry

        # Compare N to actual N and throw an error if it doesn't match
        dt_tmp <- add_N_using_filename(file_info, dt_list[[file]])
        dt_list[[file]] <- dt_tmp %>% filter(Group != "Cauchy")
    }

    dt <- rbindlist(dt_list, use.names=TRUE)
    dt <- dt %>% filter((!is.na(Pvalue)) & (!is.na(Pvalue_SKAT)) & (!is.na(Pvalue_Burden)))

    # Nudge P-values that are close to exactly 1 (for Stouffer).
    dt <- dt %>% mutate(
        Pvalue = ifelse(Pvalue > 0.99, 0.99, Pvalue),
        Pvalue_SKAT = ifelse(Pvalue_SKAT > 0.99, 0.99, Pvalue_SKAT),
        Pvalue_Burden = ifelse(Pvalue_Burden > 0.99, 0.99, Pvalue_Burden)
    )

    dt_meta <- list()
    for (test in tests)
    {
        dt_meta[[test]] <- list()
        Pvalue_col <- ifelse(test == "SKAT-O", "Pvalue", paste0("Pvalue_", test))
        
        # Weighted Fisher's meta-analysis of p-values
        dt_meta[[test]][["weighted Fisher"]] <- run_weighted_fisher(
            dt %>% group_by(Region, Group, max_MAF),
            "N_eff", Pvalue_col, "Pvalue",
            two_tail = ifelse(test == "Burden", TRUE, FALSE),
            input_beta = ifelse(test == "Burden", "BETA_Burden", NULL)) %>% 
        mutate(Stat = NA, type="weighted Fisher")

        # Stouffer's Z - Make sure P-values match, Stat= weighted_Z_Burden_Stouffer
        dt_meta[[test]][["Stouffer"]] <- run_stouffer(dt %>% group_by(Region, Group, max_MAF),
            "N_eff", "Stat", Pvalue_col, "Pvalue",
            two_tail = ifelse(test == "Burden", TRUE, FALSE),
            input_beta = ifelse(test == "Burden", "BETA_Burden", NULL)) %>% 
        mutate(type="Stouffer")
        dt_meta[[test]] <- rbindlist(dt_meta[[test]], use.names=TRUE) %>% mutate(class=test)
    }

    dt_meta <- rbindlist(dt_meta)
    print(dt_meta)

    dt_cauchy <- list()
    dt_n_eff <- unique(dt %>% select(Region, dataset, ancestry, N_eff))
    setkeyv(dt_n_eff, c("Region", "dataset", "ancestry", "N_eff"))
    for (test in tests) {
        cat(paste0(test, "..."))
        Pvalue_col <- ifelse(test == "SKAT-O", "Pvalue", paste0("Pvalue_", test))        
        dt_cauchy[[test]] <- run_cauchy(dt %>% group_by(Region, dataset, ancestry),
            "N_eff", "Stat", Pvalue_col, "Pvalue") %>% mutate(type="Cauchy")
        dt_cauchy[[test]] <- merge(dt_cauchy[[test]], dt_n_eff)
    }
    cat("\n")

    dt_meta_cauchy <- list()
    for (test in tests)
    {
        dt_meta_cauchy[[test]] <- list()
        # Weighted Fisher's meta-analysis of p-values
        dt_meta_cauchy[[test]][["weighted Fisher"]] <- run_weighted_fisher(
            dt_cauchy[[test]] %>% group_by(Region),
            "N_eff", "Pvalue", "Pvalue") %>% mutate(Stat = NA, type="weighted Fisher")
        # Stouffer's Z - Make sure P-values match, Stat = weighted_Z_Burden_Stouffer
        dt_meta_cauchy[[test]][["Stouffer"]] <- run_stouffer(dt_cauchy[[test]] %>% group_by(Region),
            "N_eff", "Stat", "Pvalue", "Pvalue") %>% mutate(type="Stouffer")
        dt_meta_cauchy[[test]] <- rbindlist(dt_meta_cauchy[[test]], use.names=TRUE) %>% mutate(class=test)
    }
    dt_meta_cauchy <- rbindlist(dt_meta_cauchy)
    dt_meta <- rbind(dt_meta, dt_meta_cauchy %>% mutate(Group="Cauchy", max_MAF="Cauchy"))
    fwrite(dt_meta, file=ifelse(grepl(".tsv.gz$", args$out), args$out, paste0(args$out, ".tsv.gz")), sep='\t')

    # Finally, carry out heterogeneity test
    dt_het <- run_heterogeneity(
        dt %>% group_by(Region, Group, max_MAF),
        "N_eff", "BETA_Burden", "BETA_Burden_meta")

}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--file_paths", default=NULL, required=TRUE,
    help="Comma separated file paths for files for meta-analysis")
parser$add_argument("--case_control_threshold", default=100, required=FALSE,
    help=paste0("Case/control count threshold for inclusion of the data into the meta-analysis [default=100]"))
parser$add_argument("--out", default="meta_analysis", required=FALSE,
    help="Output filepath")
args <- parser$parse_args()

main(args)
