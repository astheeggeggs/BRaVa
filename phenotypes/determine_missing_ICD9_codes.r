library(data.table)
library(dplyr)
library(googlesheets4)
library(readxl)

source("BRaVa_phenotypes_utils.r")

dt <- munge_BRaVa_ICD_proposals()

# ICD10 codes to map
dt_ICD9_case_required <- dt %>% filter((ICD10_case_include != "") & (ICD9_case_include == ""))
dt_ICD9_control_required <- dt %>% filter((ICD10_control_exclude != "") & (ICD9_control_exclude == ""))

to_match <- list(case = dt_ICD9_case_required, control = dt_ICD9_control_required)

for (entry in to_match) {
    for (i in 1:nrow(entry)) {
        strings_to_map <- strsplit(gsub("^ (.*) $", "\\1", entry$ICD10_case_include[i]), split=" \\| ")[[1]]
        # If the ICD10 code is 2 digit, add a * as the dictionary matches to 3 digit
        ICD_digit <- nchar(strings_to_map)
        if (any(ICD_digit == 3)) {
            strings_to_map[which(ICD_digit == 3)] <- paste0(strings_to_map[which(ICD_digit == 3)], ".*")
        }
        strings_to_map <- paste0(strings_to_map, collapse="|")
        ICD9_case <- extract_terms(strings_to_map, grep=TRUE)
        cat(paste0(entry$phenotype[i], " (phenotype ID: ", entry$phenotypeID[i], ")\n"))
        cat(paste0(paste0(strings_to_map, collapse=", "), "\n"))
        cat(paste0(paste0("-> ", paste0(ICD9_case, collapse=", ")), "\n"))
    }
}
