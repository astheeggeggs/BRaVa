# install.packages("googlesheets4")
library(data.table)
library(dplyr)
library(ggplot2)
library(googlesheets4)
library(stringr)

source("BRaVa_phenotypes_utils.r")

# Extract continuous traits
extract_continuous_trait_counts()

dt_query_ICD <- munge_BRaVa_ICD_proposals()
fwrite(dt_query_ICD, file='data/BRaVa_ICD_proposals.tsv', sep='\t', quote=TRUE)
dt_query_ICD <- fread('data/BRaVa_ICD_proposals.tsv')
dt_data <- extract_ID_and_ICD_UKB('data/ukb10844_ukb50009_updateddiagnoses_14012022.csv.gz')
dt_binary_list <- extract_case_status(dt_data, dt_query_ICD)
dt_binary_ICD <- dt_binary_list$dt_data
dt_query_ICD <- dt_binary_list$dt_query
fwrite(dt_query_ICD, file='data/BRaVa_ICD_proposals_grep.tsv', sep='\t', quote=TRUE)
setkey(dt_binary_ICD, "eid")

# Merge with 1000G labels - this file is created using 05_estimate_superpopulation.r in the QC folder.
dt_classify <- fread("/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv")
dt_classify[, eid:=sample.ID]
dt_classify[, sample.ID:=NULL]
setkey(dt_classify, "eid")

dt_binary_classified <- merge(dt_binary_ICD, dt_classify)

dt_query_OPCS4 <- munge_BRaVa_OPCS4_proposals()
fwrite(dt_query_OPCS4, file='data/BRaVa_OPCS4_proposals.tsv', sep='\t', quote=TRUE)
dt_query_OPCS4 <- fread('data/BRaVa_OPCS4_proposals.tsv')
dt_data <- extract_ID_and_OPCS4_UKB('data/ukb10844_ukb50009_updateddiagnoses_14012022.csv.gz')

dt_binary_procedures <- extract_procedure_case_status(dt_data, dt_query_OPCS4)
ICD_and_procedures <- intersect(names(dt_binary_procedures)[-1], names(dt_binary_ICD)[-1])
if (length(ICD_and_procedures) >= 1) {
	where <- which(names(dt_binary_procedures) %in% ICD_and_procedures)
	names(dt_binary_procedures)[where] <- paste0(names(dt_binary_procedures)[where], " procedures")
}
setkey(dt_binary_procedures, "eid")

dt_binary_classified <- merge(dt_binary_classified, dt_binary_procedures)

for (i in 1:length(ICD_and_procedures)) {
	dt_binary_classified[[paste(ICD_and_procedures[i], "ICD")]] <- dt_binary_classified[[ICD_and_procedures[i]]]
	dt_binary_classified[[paste(ICD_and_procedures[i])]][which(dt_binary_classified[[paste(ICD_and_procedures[i], "procedures")]] == 1)] <- 1
}

fwrite(dt_binary_classified, file="/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/BRaVa_phenotypes_with_superpopulation_labels_updated.tsv", sep="\t")

# Split by 1000G label and count
dt_counts <- dt_binary_classified %>% group_by(classification_strict) %>% summarise(across(c(union(dt_query_ICD$phenotype, dt_query_OPCS4$phenotype), paste(ICD_and_procedures, c("procedures", "ICD"))), sum, na.rm=TRUE))
dt_counts_t <- data.table::transpose(dt_counts, keep.names="phenotype", make.names="classification_strict")

# Finally, comma separate and combine the ICD9 and ICD10 codes together for inclusion on ICD_phecode tab of BRaVa_Nominate_Phenotypes spreadsheet.
fwrite(dt_counts_t, file="data/output/UKBB_case_counts_updated.tsv", sep="\t")

dt_counts <- dt_binary_classified %>% summarise(across(c(union(dt_query_ICD$phenotype, dt_query_OPCS4$phenotype), paste(ICD_and_procedures, c("procedures", "ICD"))), sum, na.rm=TRUE))
dt_counts_t <- data.table::transpose(dt_counts, keep.names="phenotype")
names(dt_counts_t)[2] <- "count"
fwrite(dt_counts_t, file="data/output/UKBB_case_counts_total_updated.tsv", sep="\t")

# To fill in the google sheets table
dt_counts <- merge(fread("data/output/UKBB_case_counts_updated.tsv"), fread("data/output/UKBB_case_counts_total_updated.tsv"))

# Merge in the age, sex, and age*sex, age^2*sex covariates for the continuous and case control traits.
dt_cov <- extract_covariates('data/ukb10844_ukb50009_updateddiagnoses_14012022.csv.gz')
dt_binary_classified <- fread("/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/BRaVa_phenotypes_with_superpopulation_labels_updated.tsv")
dt_cts_classified <- fread("/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/BRaVa_cts_phenotypes_with_superpopulation_labels_updated.tsv")

dt_cov$eid <- as.character(dt_cov$eid)
dt_binary_classified$eid <- as.character(dt_binary_classified$eid)
dt_cts_classified$eid <- as.character(dt_cts_classified$eid)

setkey(dt_cov, "eid")
cols <- c(
	"ICD9_string",
	"ICD10_string",
	"classification_loose",
	"super.population",
	"population",
	"OPCS4_string")
dt_binary_classified <- dt_binary_classified[, -cols, with=FALSE]
setkey(dt_binary_classified, "eid")

cols <- c(paste0("PC", seq(1,10)),
	"super.population", "population",
	"classification_strict", "classification_loose"
	)
dt_cts_classified <- dt_cts_classified[, -cols, with=FALSE]
setkey(dt_cts_classified, "eid")

dt <- merge(dt_cov, dt_binary_classified)
dt <- merge(dt, dt_cts_classified)

# Clean up names
names(dt) <- gsub(" [ ]+", " ", names(dt))
names(dt) <- gsub(" / ", "/", names(dt))
names(dt) <- gsub("([a-z])(\\()", "\\1 \\2", names(dt))
names(dt) <- gsub(" \\([A-Z,\\/a-z]+\\)", "", names(dt))
names(dt) <- gsub(",", "", names(dt))
names(dt) <- gsub("-", " ", names(dt))
names(dt) <- gsub("[\\(,\\)]", "", names(dt))
names(dt) <- gsub(" ", "_", names(dt))

# Add in IID
dt <- dt %>% rename(IID = eid)
dt <- dt %>% mutate(age=as.integer(age), sex=as.integer(sex)) %>% mutate(age_sex = age*sex, age2 = age^2, age2_sex=age^2*sex)

fwrite(dt, file="/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/BRaVa_phenotypes_with_superpopulation_labels_updated_combined.tsv", sep='\t')
# Male only phenotype information
fwrite(dt %>% filter(sex == 1) %>% select(-c("age_sex", "age2_sex")), file="/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/BRaVa_phenotypes_with_superpopulation_labels_updated_combined_M.tsv")
# Female only phenotype information
fwrite(dt %>% filter(sex == 0) %>% select(-c("age_sex", "age2_sex")), file="/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/BRaVa_phenotypes_with_superpopulation_labels_updated_combined_F.tsv")
