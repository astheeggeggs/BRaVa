library(data.table)
library(dplyr)

# Find the METAL files to run, and loop over the meta-analyses
METAL_file_location <- "~/Repositories/BRaVa_curation/data/meta_analysis/meta_results/variant"
METAL_files <- grep("^.*/METAL_.*.txt$", dir(METAL_file_location, full.names=TRUE), value=TRUE)

METAL_location <- "~/Repositories/METAL/build/bin/metal"
# Use intall_metal.sh script is METAL is not installed

for (file in METAL_files) {
	cat(paste0("carrying out meta-analysis using METAL_file: ", file, "\n"))
	system(paste(METAL_location, file))
}

cat("completed meta-analysis of all traits.\n")
