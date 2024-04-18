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

# Expected annotation names
correct_names <- c("pLoF",
	"damaging_missense_or_protein_altering",
	"other_missense_or_protein_altering",
	"synonymous",
	"pLoF;damaging_missense_or_protein_altering",
	"pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
	"Cauchy")

requested_max_MAFs <- c(1e-4, 1e-3, 1e-2)

# Expected columns of gene results files
default_gene_result_columns <- c(
	"Region", "Group", "max_MAF", "Pvalue",
	"Pvalue_Burden", "Pvalue_SKAT","BETA_Burden", "SE_Burden",
	"MAC", "MAC_case", "MAC_control", "Number_rare", "Number_ultra_rare")

minimal_gene_result_columns <- c(
    "Region", "Group", "max_MAF", "Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
    "BETA_Burden", "SE_Burden")

renaming_plot_group_list <- list(
    damaging_missense_or_protein_altering = "Damaging missense or PA",
    other_missense_or_protein_altering = "Other missense or PA",
    synonymous = "Synonymous",
    pLoF = "pLoF"
)

file_check_information <- list(
	datasets = c(
		"all-of-us",
		"alspac",
		"biome",
		"bbj",
		"ckb",
		"ccpm",
		"decode",
		"egcut",
		"dan-rav",
		"genes-and-health",
		"gel",
		"mgbb",
		"pmbb",
		"qatar-genomes",
		"uk-biobank",
		"viking-genes"
	),
	phenotypes = c(
		"AAA", "AcApp","AcuLymLeuk", "Adenomy", "AMD", "ALamy", "AUD","AloAre",
		"AnoNer", "AoSten", "Asth", "AtopDis", "AFib", "ADHD", "ASD", "BCLL",
		"BenCervUterNeo", "BenIntNeo", "BenNodGoit", "BladCanc", "BrainCNSCanc",
		"BreastCanc", "BrugSynd", "BuliNer", "BullPemph", "CarShock", "HCM","CRVO",
		"CervCanc", "CML", "COPD", "CRF", "CoffSirSynd", "ColonRectCanc", "CAD",
		"CCANS", "EatDis", "Endocar", "Endometr", "EsophCanc", "EssThrom", "EFRMB",
		"FSP", "FemInf", "FemInfAC", "FolLymph", "Gout", "GravesDis", "HemoChromo",
		"HF", "HepCarcin","HTN", "HHD", "HypoThyr", "HypoThyrSec", "IPF", "ITP",
		"IBD", "IFHern", "ILDSarc", "IodDef", "KabSynd", "KidCanc", "KleefSynd",
		"LaryxCanc", "Leuk", "LiverCanc", "LiverFibCirr", "LongQTSynd",
		"LymphThyrit", "MalInf", "MatHem", "MatHypDis", "MODYDiab", "MultiMyel",
		"MS", "MECS", "Myocard", "Narco1", "NonFuncPitAd", "NHL", "NonPapTCCBlad",
		"NonRheuValv", "OUD", "OCD", "OvCanc", "Pancreat", "ParkDis", "PeptUlcer",
		"PAD", "PlacInsuf", "PCOS", "PolycythVera", "Preeclamps", "PregLoss",
		"POAG", "PrimSjoSynd", "Prolactinom", "Psori", "RheumHeaDis", "RheumArth",
		"RomWardSynd", "Sarcoid", "SebDerm", "SpinaBifAp", "StomCanc", "Stroke",
		"SLE", "TAAD", "ThyroCanc", "T2Diab", "Urolith", "UterCanc", "VaricVeins",
		"VTE", "ALT", "AlcCons", "AST", "BMI", "CRP", "CACS", "CK", "HDLC",
		"Height", "LDLC", "TChol", "TG", "WHRBMI", "LVH", "Append", "HipRep",
		"CogAbil", "EduAtt", "PsySymp", "SchGrades", "SCDCAT"
	),
	sexes = c(
		"ALL",
		"M",
		"F"
	),
	ancestries = c(
		"AFR",
		"AMR",
		"EAS",
		"EUR",
		"MID",
		"SAS"
	),
	types = c(
		"gene",
		"variant"
	)
)

search_for_files <- function(file)
{
    if (file.exists(paste0("gene/", file, ".txt.gz")) &
        file.exists(paste0("variant/", file, ".txt.singleAssoc.txt.gz"))) {
        gene_file <- paste0("gene/", file, ".txt.gz")
        variant_file <- paste0("variant/", file, ".txt.singleAssoc.txt.gz")
        gz <- TRUE
    } else if (file.exists(paste0("gene/", file, ".txt")) &
               file.exists(paste0("variant/", file, ".txt.singleAssoc.txt"))) {
        gene_file <- paste0("gene/", file, ".txt")
        variant_file <- paste0("variant/", file, ".txt.singleAssoc.txt")
        gz <- FALSE
    } else if (file.exists(paste0("gene/", file, "_F.txt.gz")) &
               file.exists(paste0("variant/", file, "_F.txt.singleAssoc.txt.gz"))) {
        gene_file <- paste0("gene/", file, "_F.txt.gz")
        variant_file <- paste0("variant/", file, "_F.txt.singleAssoc.txt.gz")
        gz <- TRUE
    } else if (file.exists(paste0("gene/", file, "_F.txt")) &
               file.exists(paste0("variant/", file, "_F.txt.singleAssoc.txt"))) {
        gene_file <- paste0("gene/", file, "_F.txt")
        variant_file <- paste0("variant/", file, "_F.txt.singleAssoc.txt")
        gz <- FALSE
    } else {
        return(NULL)
    }
    return(list(gene_file=gene_file, variant_file=variant_file, gz=gz))
}

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

extract_file_info <- function(filename)
{
	gz <- ifelse(grepl(".gz$", filename), TRUE, FALSE)
	filename <- gsub(".gz$", "", filename)
	filename <- gsub("cleaned.", "", filename)
	file_info <- as.list(strsplit(filename, split="\\.")[[1]])

	if (!(length(file_info) %in% c(12, 13))) {
		print(file_info)
		stop("Incorrect file naming convention, please check filenames")
	}

	if (length(file_info) == 13) {
		binary <- TRUE
		names(file_info) <- c(
			"dataset",
			"last_name",
			"analysis_name",
			"phenotype",
			"freeze_number",
			"sex",
			"ancestry",
			"n_cases",
			"n_controls",
			"software",
			"type",
			"date",
			"split")
	} else {
		binary <- FALSE
		names(file_info) <- c(
			"dataset",
			"last_name",
			"analysis_name",
			"phenotype",
			"freeze_number",
			"sex",
			"ancestry",
			"n",
			"software",
			"type",
			"date",
			"split")
	}

	file_info$gz <- gz
	file_info$binary <- binary
	return(file_info)
}

add_N_using_filename <- function(file_info, dt)
{
    if (!file_info$binary) {
        dt$N_eff <- file_info$n
    } else {
        N_case <- file_info$n_cases
        N_control <- file_info$n_controls
        N_eff <- (4 * N_case * N_control) / (N_case + N_control)
        dt$N_eff <- N_eff
        dt$N_case <- N_case
        dt$N_control <- N_control
    }
    return(dt)
}

add_N <- function(file_info, dt_gene)
{
    # Read in the variant file to extract the sample size
    dt_variant <- fread(cmd = ifelse(file_info$gz,
        paste0("gzcat ", file_info$variant_file),
        paste0("cat ", file_info$variant_file)))

    if ("N" %in% colnames(dt_variant)) {
        N_eff <- mean(dt_variant$N)
        dt_gene$N_eff <- N_eff
        binary <- FALSE
    } else {
        N_case <- mean(dt_variant$N_case)
        N_control <- mean(dt_variant$N_ctrl)
        N_eff <- (4 * N_case * N_control) / (N_case + N_control)
        dt_gene$N_eff <- N_eff
        dt_gene$N_case <- N_case
        dt_gene$N_control <- N_control
        binary <- TRUE
    }
    return(list(dt_gene = dt_gene, binary=binary))
}

cauchy_combination <- function(p_values, weights=NULL)
{
	is.zero <- sum(p_values == 0) >= 1
	is.one <- sum(p_values > (1 - 1e-14)) >= 1

	if (is.zero) {
		return(0)
	}

	if (is.one) {
		p <- min(p_values) * length(p_values)
		if (p > 1) {
			return(-Inf)
		} else {
			return(qcauchy(p, lower.tail=FALSE))
		}
		return(min(1, (min(p_values)) * (length(p_values))))
	}

	if (is.null(weights)) {
		weights <- rep(1 / length(p_values), length(p_values))
	} else if (length(weights) != length(p_values)) {
		stop("The length of weights should be the same as that of the p-values!")
	} else if (sum(weights < 0) > 0) {
		stop("All the weights must be positive!")
	} else {
		weights <- weights / sum(weights)
	}

	is_small <- (p_values < 1e-16)
	if (sum(is_small) == 0) {
		cct_stat <- sum(weights * tan((0.5 - p_values) * pi))
	} else {
		cct_stat <- sum((weights[is_small] / p_values[is_small]) / pi)
		cct_stat <- cct_stat + 
			sum(weights[!is_small] * tan((0.5 - p_values[!is_small]) * pi))
	}

	return(cct_stat)
}

weighted_fisher <- function(p_values, weights=NULL, two_tail=FALSE, input_beta=NULL)
{
	if (is.null(weights)) {
		weights <- rep(1, length(p_values))
	}

	idx.na <- which(is.na(p_values))
	
	if (length(idx.na) > 0) {
		p_values <- p_values[-idx.na]
		weights <- weights[-idx.na]
		if (two_tail) {
			input_beta <- input_beta[-idx.na]
		}
	}

	NP <- length(p_values)
	NS <- length(weights)
	if (NP != NS) { stop("The length of p and weights vector must be identical.") }

	N <- NS
	Ntotal <- sum(weights)
	ratio <- weights / Ntotal
	Ns <- N * ratio
	G <- c()

	if (!two_tail) {
		for (i in 1:length(p_values)) {
		  G <- append(G, qgamma(p = p_values[i], shape = Ns[i], scale = 2, lower.tail=FALSE))
		}
		Gsum <- sum(G)
		resultP <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = FALSE)
	} else {
		p1 <- p2 <- p_values
		idx_pos <- which(input_beta > 0)
		idx_neg <- which(input_beta < 0)
		
		# Positive direction
		G <- c()
		p1[idx_pos] <- p_values[idx_pos] / 2
		p1[idx_neg] <- 1 - p_values[idx_neg] / 2

		for (i in 1:length(p1)) {
		  G <- append(G, qgamma(p = p1[i], shape = Ns[i], scale = 2, lower.tail = FALSE))
		}
		Gsum <- sum(G)
		resultP1 <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = FALSE)
		
		# Negative direction
		G <- c()
		p2[idx_pos] <- 1 - p_values[idx_pos] / 2
		p2[idx_neg] <- p_values[idx_neg] / 2

		for (i in 1:length(p2)) {
		  G <- append(G, qgamma(p = p2[i], shape = Ns[i], scale = 2, lower.tail = FALSE))
		}
		Gsum <- sum(G)
		resultP2 <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = FALSE)
		resultP <- 2 * min(resultP1, resultP2)
		if (resultP > 1.0) {
			resultP <- 1.0
		}
	}
	return(min(1, resultP))
}

run_weighted_fisher <- function(
	grouped_dt, n_eff_name, input_pvalues, output_meta_pvalue,
	two_tail = FALSE, input_beta = NULL)
{
	if (two_tail) {
		result <- grouped_dt %>% 
		summarise(
			"{output_meta_pvalue}" := weighted_fisher(
				.data[[input_pvalues]],
				weights = if (is.null(n_eff_name)) { NULL } else { .data[[n_eff_name]] },
				two_tail = two_tail, input_beta=.data[[input_beta]]
			)
		)
	} else {
		result <- grouped_dt %>% 
		summarise(
			"{output_meta_pvalue}" := weighted_fisher(
				.data[[input_pvalues]],
				weights = if (is.null(n_eff_name)) { NULL } else { .data[[n_eff_name]] }
			)
		)
	}
	return(result)
}

create_gene_data_table <- function(group_files_lines)
{
	group_files_variant_line <- group_files_lines[1]
	group_files_group_line <- group_files_lines[2]
	to_munge_variant <- strsplit(group_files_variant_line, split=" ")[[1]]
	n <- length(to_munge_variant)
	to_munge_group <- strsplit(group_files_group_line, split=" ")[[1]]
	if (n != length(to_munge_group)) {
		stop("Error: number of variants in the region does not match the number of groups in the region")
	}
	return(data.table(Region=to_munge_variant[1],
		MarkerID=to_munge_variant[3:n],
		Group=to_munge_group[3:n]))
}

extract_group_file_information <- function(regexp_group_files) {
	# Extact the directory portion and the regular expression within the folder
	regexp_group_files <- strsplit(regexp_group_files, "/")[[1]]
	if (length(regexp_group_files) > 1) {
		directory <- paste(regexp_group_files[1:(length(regexp_group_files)-1)], collapse="/")
		regexp_group_files <- regexp_group_files[length(regexp_group_files)]
	} else {
		directory <- "."
	}
	group_files <- lapply(dir(directory, pattern=regexp_group_files), fread, sep="@", header=FALSE)
	group_files <- rbindlist(group_files)
	print(group_files)
	gene_data_table_list <- list()
	total_lines <- nrow(group_files)
	i <- 0
	for (i in 1:(nrow(group_files)/2)) {
		print(i)
		gene_data_table_list[[i]] <- create_gene_data_table(group_files$V1[(2*i-1):(2*i)])
	}
	gene_data_table <- rbindlist(gene_data_table_list)
	return(gene_data_table)
}

combine_gene_and_variant_information <- function(variant_file, gene_data_table)
{
	# Merge with the variant files to obtain the variants that have a MAC>10 in the analysed dataset
	dt_variant <- fread(variant_file)
	# Deal with the ultra-rare variants
	setkey(dt_variant, "MarkerID")
	setkey(gene_data_table, "MarkerID")

	dt <- merge(dt_variant, gene_data_table, all=TRUE)
	dt[, Region:=ifelse(CHR=="UR", gsub(":.*", "", MarkerID), Region)]
	dt[, Group:=ifelse(CHR=="UR", gsub(".*:(.*):.*", "\\1", MarkerID), Group)]
	dt[, Tstat := ifelse(AF_Allele2 > 0.5, -Tstat, Tstat)]
	dt[, MAF := pmin(AF_Allele2, 1-AF_Allele2)]
	return(dt)
}

extract_combined_annotations <- function(gene_file) {
	dt_gene <- fread("data/meta_analysis/gcloud/uk-biobank.palmer.PRELIMINARY.Coronary_artery_disease.JULY23Freeze.ALL.EUR.19915.382460.SAIGE.gene.20240110.cleaned.txt.gz")
	combined_annotations <- grep(";", unique(dt_gene$Group), value=TRUE)
	return(combined_annotations)
}

extract_beta_burden <- function(dt, combined_annotations, max_MAF)
{
	dt <- dt %>% filter(!is.na(Region))
	m <- max_MAF
	dt_common <- dt %>% 
		filter(CHR != "UR", MAF < max_MAF) %>% 
		group_by(Region, Group) %>% 
		summarise(
			numerator_common = sum(Tstat),
			denominator_common = sum(var))
	dt_rare <- dt %>% filter(CHR == "UR") %>% 
		mutate(max_MAF = as.numeric(gsub(".*:", "", MarkerID))) %>% 
		filter(max_MAF == m) %>% group_by(Region, Group) %>% 
		summarise(
			numerator_rare = sum(Tstat),
			denominator_rare = sum(var))
	dt_common <- data.table(dt_common)
	
	# Finally, if there are any combined annotations, then
	# include those

	# The below just needs to be done for dt_common, not for dt_rare.
	dt_common_combined_list <- list()
	for (str_comb_annot in combined_annotations) {
		combined_annotation_vector <- strsplit(str_comb_annot, split=";")[[1]]
		dt_common_combined_list[[str_comb_annot]] <- dt_common %>% group_by(Region) %>% 
			filter(Group %in% combined_annotation_vector) %>%
			summarise(
				numerator_common = sum(numerator_common),
				denominator_common = sum(denominator_common)
				) %>% mutate(Group = str_comb_annot)
	}
	rbind(rbindlist(dt_common_combined_list), dt_common)

	setkeyv(dt_common, c("Region", "Group"))
	dt_rare <- data.table(dt_rare)
	setkeyv(dt_rare, c("Region", "Group"))
	dt <- merge(dt_common, dt_rare, all=TRUE)
	
	# Numerator
	dt[, numerator_common := ifelse(is.na(numerator_common), 0, numerator_common)]
	dt[, numerator_rare := ifelse(is.na(numerator_rare), 0, numerator_rare)]

	# Denominator
	dt[, denominator_common := ifelse(is.na(denominator_common), 0, denominator_common)]
	dt[, denominator_rare := ifelse(is.na(denominator_rare), 0, denominator_rare)]
	return(dt)
}

extract_beta_burden_across_MAFs <- function(dt, max_MAFs, combined_annotations) {
	# Loop over max MAFs
	dt_list <- list()
	# unique(dt_gene$max_MAF[!is.na(dt_gene$max_MAF)])
	for (m in max_MAFs) {
		dt_list[[as.character(m)]] <- extract_beta_burden(dt, combined_annotations, m)
		dt_list[[as.character(m)]]$max_MAF <- m
	}
	dt_out <- rbindlist(dt_list) %>% group_by(Region, Group, max_MAF)
}

run_heterogeneity <- function(
	grouped_dt, n_eff_name, input_beta, output_meta_beta)
{
	grouped_dt <- grouped_dt %>%
	mutate(
		weights = sqrt(.data[[n_eff_name]]),
		beta = .data[[input_beta]]
	)
	summary_dt <- grouped_dt %>% 
	summarise(
		sum_weights = sum(weights),
		sum_betas = sum(weights * beta)
		) %>% mutate("{output_meta_beta}" := sum_betas/sum_weights)
	return(
		merge(grouped_dt, summary_dt) %>% 
		mutate(
			deviation = weights * (beta - .data[[output_meta_beta]])^2
		) %>% group_by(across(as.character(groups(grouped_dt)))) %>%
		summarise(sum_deviation = sum(deviation), n_studies=n()) %>% 
		filter(n_studies > 1) %>% mutate(p_het = pchisq(sum_deviation, n_studies-1, lower.tail=FALSE))
	)
}

run_fisher <- function(
	grouped_dt, chi2_stat_name, input_pvalues, output_meta_pvalue)
{
	return(
		grouped_dt %>% 
		summarise(df = 2*n(),
			"{chi2_stat_name}" := -2*sum(log(.data[[input_pvalues]]))) %>%
		mutate("{output_meta_pvalue}" := pchisq(.data[[chi2_stat_name]],
			df=df, lower.tail=FALSE)) %>% select(-df)
	)
}

run_cauchy <- function(
	grouped_dt, n_eff_name, Cauchy_stat_name, input_pvalues, output_meta_pvalue) 
{
	return(
		grouped_dt %>% 
		summarise(
			"{Cauchy_stat_name}" := cauchy_combination(
				.data[[input_pvalues]],
				weights=sqrt(.data[[n_eff_name]])
			),
			number_of_pvals := n()
		) %>% mutate(
		  	"{output_meta_pvalue}" := ifelse(
				.data[[Cauchy_stat_name]] > 1e+15,
				(1 / .data[[Cauchy_stat_name]]) / pi,
				pcauchy(.data[[Cauchy_stat_name]], lower.tail=FALSE)
			)
		) %>% mutate(
			"{output_meta_pvalue}" := ifelse(
				.data[[output_meta_pvalue]] > (1 - 1e-10),
				(1 - 1/number_of_pvals), .data[[output_meta_pvalue]]
			)
		)
	)
}

run_inv_var <- function(
	grouped_dt, beta_meta_name, se_meta_name,
	input_beta, input_SE, output_meta_pvalue
) {
	return(
		grouped_dt %>% 
		mutate(weight = 1/(.data[[input_SE]]**2)) %>%
		mutate(effs_inv_var = .data[[input_beta]] * weight) %>%
		summarise(
			"{beta_meta_name}" := sum(effs_inv_var) / sum(weight),
			"{se_meta_name}" := sqrt(1/sum(weight)),
			"{output_meta_pvalue}" := 2 * dnorm(
				abs(.data[[beta_meta_name]] / .data[[se_meta_name]]))
		)
	)
}

run_stouffer <- function(
	grouped_dt, n_eff_name, weighted_Z_name,
	input_pvalues, output_meta_pvalue,
	two_tail = FALSE, input_beta = NULL
) {
	if (two_tail) {
		grouped_dt <- grouped_dt %>% 
			mutate("{input_pvalues}" := .data[[input_pvalues]]/2)
	} else {
		input_beta <- "beta_dummy"
		grouped_dt <- grouped_dt %>% mutate("{input_beta}" := 1)
	}
	
	result <- grouped_dt %>%
		mutate(
			weighted_Z_numerator = (
				sqrt(.data[[n_eff_name]]) * 
				(-qnorm(.data[[input_pvalues]])) * 
				sign(.data[[input_beta]])
			)
		)

	if (two_tail) {
		result <- result %>%
		summarise(
			"{weighted_Z_name}" := sum(weighted_Z_numerator) / 
				sqrt(sum(.data[[n_eff_name]])),
			"{output_meta_pvalue}" := 2 * pnorm(abs(.data[[weighted_Z_name]]), lower.tail=FALSE)
		)
	} else {
		result <- result %>%
		summarise(
			"{weighted_Z_name}" := sum(weighted_Z_numerator) / 
				sqrt(sum(.data[[n_eff_name]])),
			"{output_meta_pvalue}" := pnorm(.data[[weighted_Z_name]], lower.tail=FALSE)
		)
	}
	return(result)
}