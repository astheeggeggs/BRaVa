library(data.table)
library(dplyr)

# Find the METAL files to run, and loop over the meta-analyses
METAL_file_location <- ""
METAL_files <- dir(METAL_file_location)

METAL_location <- "~/Repositories/METAL/bin/metal"
# Use intall_metal.sh script is METAL is not installed

for (file in METAL_files) {
	cat(paste0("carrying out meta-analysis using METAL_file: ", file))
	system(paste(METAL_location, file))
}

cat("completed meta-analysis of all traits.\n")
