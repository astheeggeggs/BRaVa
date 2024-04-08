#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(latex2exp)
library(argparse)

source("~/Repositories/BRaVa_curation/QC/utils/pretty_plotting.r")
source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")

main <- function(args)
{
    files <- dir(args$analysis_results_folder,
        pattern=args$analysis_results_path_regexp,
        full.names=TRUE)
    print(args)

    files <- setdiff(files, grep("gel", files, value=TRUE))

    pdf(file=args$out, width=8, height=4)
    for (file in files)
    {
        file_info <- extract_file_info(gsub(".*/(.*)", "\\1", file))
        phe_plot <- gsub("_", " ", gsub("_$", "", str_trim(gsub("[[:space:]_]+", "\\_", file_info$phenotype))))
        phe_plot <- paste(file_info$dataset, file_info$ancestry, phe_plot, sep=", ")
        phe_plot <- paste0(phe_plot, "\n(", file_info$n_cases, " cases, ", file_info$n_controls, " controls)")
        n_cases <- as.integer(file_info$n_cases)
        n_controls <- as.integer(file_info$n_controls)
        n <- n_cases + n_controls
        cat(paste0(phe_plot, "...\n"))
        dt <- fread(cmd = paste("gzcat", file))
        dt <- melt(dt, id.vars = c("Region", "MAC", "Group", "max_MAF"),
                measure.vars = c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT"),
                value.name = "Pvalue", variable.name = "class") %>% 
        mutate(Pvalue = -log10(Pvalue))
        dt[, class := as.factor(class)]
        dt$max_MAF[which(is.na(dt$max_MAF))] <- "Cauchy"
        levels(dt$class) <- list(
            `SKAT-O` = "Pvalue",
            Burden = "Pvalue_Burden",
            SKAT = "Pvalue_SKAT")
        dt[, E_MAC_case := (n_cases / n) * MAC]
        dt[, Interval := cut(E_MAC_case, breaks=c(0, 5, 10, 100, 1000, Inf),
        labels=c("<5", "[5, 10)", "[10, 100)", "[100, 1000)", ">1000"),
        include.lowest=TRUE)]

        dt_to_plot <- data.table(
            dt %>% 
            arrange(Group, max_MAF, class, Interval, desc(Pvalue)) %>%
            select(-c(Group, max_MAF, class, Interval)),
            dt %>% 
            group_by(Group, max_MAF, class, Interval) %>% 
            arrange(desc(Pvalue)) %>% 
            summarize(
                Pvalue_expected = -log10(seq(1, n())/(n() + 1)),
        	)
        )

        max_MAF_groups <- setdiff(unique(dt_to_plot$max_MAF), "Cauchy")
        groups <- setdiff(unique(dt_to_plot$Group), "Cauchy")
        for (m in max_MAF_groups) {
            for (g in groups) {
                max_MAF_plot <- ifelse(
                    grepl("e", as.character(m)),
                    gsub("([0-9\\.]+)e(-)*0*([1-9][0-9]*)", "$\\1\\\\times 10^{\\2\\3}$", as.character(m)),
                    paste0("$", as.character(m), "$"))
                variant_class_plot <- gsub("_", " ", gsub("[\\|;]", ", ", g))
                p <- create_pretty_qq_plot(
                    plot_title=phe_plot,
                    plot_subtitle=TeX(paste0(variant_class_plot, "; max MAF = ", max_MAF_plot)),
                    cex_labels=2,
                    dt_to_plot %>% filter(Group==g, max_MAF==m),
                    aes(x=Pvalue_expected, y=Pvalue, color=Interval),
                    save_figure=FALSE,
                    x_label=TeX("$-\\log_{10}(P_{expected})$"), 
                    y_label=TeX("$-\\log_{10}(P_{observed})$"),
                    key_cols=c("class", "Pvalue"),
                    include_qq_ribbon=FALSE,
                    width=170,
                    height=120,
                    by_chr=FALSE,
                    print_p=FALSE
                ) + facet_wrap(~class)
                print(p)
            }
        }
    }
    dev.off()
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--analysis_results_folder",
    default="~/Repositories/BRaVa_curation/data/meta_analysis/gcloud", required=FALSE,
    help="Folder containing the analysis results files")
parser$add_argument("--analysis_results_path_regexp", default="cleaned", required=FALSE,
    help="Regular expression of analysis results files to loop over")
parser$add_argument("--out", default="~/Repositories/BRaVa_curation/plots/meta_analysis/analysis_all_biobanks_qq_MAC.pdf", required=FALSE,
    help="Output filepath")
args <- parser$parse_args()

# Debug
# args <- list()
# args$analysis_results_folder <- "../data/meta_analysis/gcloud"
# args$analysis_results_path_regexp <- "*cleaned*"

main(args)
