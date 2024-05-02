#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(latex2exp)
library(argparse)







gcloud_bucket <- paste0("brava-meta-upload-", dataset)
# Upload the results to the relevant bucket
system(paste0("gsutil cp ", variant_folder, "/combined/* gs://", gcloud_bucket, "/"))
system(paste0("gsutil cp ", gene_folder, "/combined/* gs://", gcloud_bucket, "/"))


# The following loop is to combine together files for meta-analysis
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


# "AFib"
# "AMD"
# "Asth"
# "BenCervUterNeo"
# "BenIntNeo"
# "BMI"
# "BreastCanc"
# "CAD"
# "CervCanc"
# "ColonRectCanc"
# "COPD"
# "CRF"
# "CRP"
# "EFRMB"
# "FemInf"
# "Gout"
# "HDLC"
# "Height"
# "HF"
# "HTN"
# "IBD"
# "IFHern"
# "ILDSarc"
# "LDLC"
# "MatHem"
# "NonRheuValv"
# "PAD"
# "Pancreat"
# "PeptUlcer"
# "Psori"
# "RheumArth"
# "RheumHeaDis"   
# "Stroke"
# "T2Diab"
# "TChol"
# "TG"
# "Urolith"
# "VaricVeins"
# "VTE"
# "WHRBMI"
"AcApp"
# "AlcCons"
# "ALT"
# "AST"
"CK"
"CRVO"
# "HipRep"
