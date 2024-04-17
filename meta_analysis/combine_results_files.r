#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(latex2exp)
library(argparse)

source("~/Repositories/BRaVa_curation/QC/utils/pretty_plotting.r")
source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")

analysis_results_folder <- "~/Repositories/BRaVa_curation/data/meta_analysis/gcloud"
analysis_results_path_regexp <- "*cleaned*"
files <- dir(analysis_results_folder,
    pattern=analysis_results_path_regexp,
    full.names=TRUE)

dt_list <- list()
for (file in files) {
    file_info <- extract_file_info(gsub(".*/(.*)", "\\1", file))
    print(file_info)
    pheno <- file_info$phenotype
    dataset <- file_info$dataset
    ancestry <- file_info$ancestry
    dt_tmp <- fread(cmd = paste("gzcat", file))
    dt_tmp <- melt(dt_tmp, id.vars = c("Region", "Group", "max_MAF"),
            measure.vars = c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT"),
            value.name = "Pvalue", variable.name = "class") %>% 
    mutate(Pvalue = -log10(Pvalue))
    dt_tmp[, class := as.factor(class)]
    dt_tmp$max_MAF[which(is.na(dt_tmp$max_MAF))] <- "Cauchy"
    levels(dt_tmp$class) <- list(
        `SKAT-O` = "Pvalue",
        Burden = "Pvalue_Burden",
        SKAT = "Pvalue_SKAT")
    dt_tmp$pheno <- pheno
    dt_tmp$dataset <- dataset
    dt_tmp$ancestry <- ancestry
    dt_list[[file]] <- dt_tmp
}

dt <- rbindlist(dt_list)
