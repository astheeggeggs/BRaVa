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
    ribbon_p <- 0.95
    files <- dir(args$analysis_results_folder,
        pattern=args$analysis_results_path_regexp,
        full.names=TRUE)
    print(files)
    print(args)

    pdf(file=args$out, width=6, height=4)
    for (file in files) {
        file_info <- extract_file_info(gsub(".*/(.*)", "\\1", file))
        phe_plot <- gsub("_", " ", gsub("_$", "", str_trim(gsub("[[:space:]_]+", "\\_", file_info$phenotype))))
        phe_plot <- paste(file_info$dataset, file_info$ancestry, phe_plot, sep=", ")
        phe_plot <- paste0(phe_plot, "\n(", file_info$n_cases, " cases, ", file_info$n_controls, " controls)")
        cat(paste0(phe_plot, "...\n"))
        dt <- fread(cmd = paste("gzcat", file))
        dt <- melt(dt, id.vars = c("Region", "Group", "max_MAF"),
                measure.vars = c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT"),
                value.name = "Pvalue", variable.name = "class") %>% 
        mutate(Pvalue = -log10(Pvalue))
        dt[, class := as.factor(class)]
        dt$max_MAF[which(is.na(dt$max_MAF))] <- "Cauchy"
        levels(dt$class) <- list(
            `SKAT-O` = "Pvalue",
            Burden = "Pvalue_Burden",
            SKAT = "Pvalue_SKAT")

        dt_to_plot <- data.table(
            dt %>% 
            arrange(Group, max_MAF, class, desc(Pvalue)) %>%
            select(-c(Group, max_MAF, class)),
            dt %>% 
            group_by(Group, max_MAF, class) %>% 
            arrange(desc(Pvalue)) %>% 
            summarize(
                Pvalue_expected = -log10(seq(1, n())/(n() + 1)),
                clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n():1, shape1 = 1:n())),
                cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n():1, shape1 = 1:n()))
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
                    aes(x=Pvalue_expected, y=Pvalue, color=class),
                    save_figure=FALSE,
                    x_label=TeX("$-\\log_{10}(P_{expected})$"), 
                    y_label=TeX("$-\\log_{10}(P_{observed})$"),
                    key_cols=c("class", "Pvalue"),
                    aes_ribbon = aes(ymin=clower, ymax=cupper),
                    width=170,
                    height=120,
                    by_chr=FALSE,
                    print_p=FALSE
                )
                print(p)
            }
        }
        p <- create_pretty_qq_plot(
            plot_title=phe_plot,
            plot_subtitle="Cauchy",
            cex_labels=2,
            dt_to_plot %>% filter(Group == "Cauchy"),
            aes(x=Pvalue_expected, y=Pvalue, color=class),
            save_figure=FALSE,
            x_label=TeX("$-\\log_{10}(P_{expected})$"), 
            y_label=TeX("$-\\log_{10}(P_{observed})$"),
            key_cols=c("class", "Pvalue"),
            aes_ribbon = aes(ymin=clower, ymax=cupper),
            width=170,
            height=120,
            by_chr=FALSE,
            print_p=FALSE)
        print(p)
    }
    dev.off()
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--analysis_results_folder",
    default="~/Repositories/BRaVa_curation/data/meta_analysis/gcloud", required=FALSE,
    help="Folder containing the analysis results files")
parser$add_argument("--analysis_results_path_regexp", default="*cleaned*", required=FALSE,
    help="Regular expression of analysis results files to loop over")
parser$add_argument("--out", default="~/Repositories/BRaVa_curation/plots/meta_analysis/analysis_all_biobanks_qq.pdf", required=FALSE,
    help="Output filepath")
args <- parser$parse_args()

# Debug
# args <- list()
# args$analysis_results_folder <- "../data/meta_analysis/gcloud"
# args$analysis_results_path_regexp <- "*cleaned*"

main(args)
