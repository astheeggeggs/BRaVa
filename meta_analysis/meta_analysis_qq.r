#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(latex2exp)
library(argparse)

source("~/Repositories/BRaVa_curation/QC/utils/pretty_plotting.r")

main <- function(args)
{
    ribbon_p <- 0.95
    files <- dir(args$meta_analysis_results_folder, full.names=TRUE)
    files_to_include <- grep("meta", dir(args$meta_analysis_results_folder))
    files <- files[files_to_include]
    print(files)
    print(args)

    pdf(file=args$out, width=8, height=4)
    for (file in files) {
        phe <- gsub(".*/(.*)_meta.*", "\\1", file)
        phe_plot <- gsub("_", " ", gsub("_$", "", str_trim(gsub("[[:space:]_]+", "\\_", phe))))
        cat(paste0(phe_plot, "...\n"))
        dt_meta <- fread(cmd = paste("gzcat", file)) %>% mutate(Pvalue = -log10(Pvalue))
        if (!is.null(args$type)) {
            dt_meta <- dt_meta %>% filter(type == args$type)
        }
        dt_meta_to_plot <- data.table(
            dt_meta %>% 
            arrange(class, type, Group, max_MAF, desc(Pvalue)) %>%
            select(-c(class, type, Group, max_MAF)),
            dt_meta %>% 
            group_by(class, type, Group, max_MAF) %>% 
            arrange(desc(Pvalue)) %>% 
            summarize(
                Pvalue_expected = -log10(seq(1, n())/(n() + 1)),
                clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n():1, shape1 = 1:n())),
                cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n():1, shape1 = 1:n()))
        	)
        )

        max_MAF_groups <- setdiff(unique(dt_meta_to_plot$max_MAF), "Cauchy")
        groups <- setdiff(unique(dt_meta_to_plot$Group), "Cauchy")
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
                    dt_meta_to_plot %>% filter(Group==g, max_MAF==m),
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
                ) + facet_wrap(~type)
                print(p)
            }
        }

        p <- create_pretty_qq_plot(
            plot_title=phe_plot,
            plot_subtitle="Cauchy",
            cex_labels=2,
            dt_meta_to_plot %>% filter(Group == "Cauchy"),
            aes(x=Pvalue_expected, y=Pvalue, color=class),
            save_figure=FALSE,
            x_label=TeX("$-\\log_{10}(P_{expected})$"), 
            y_label=TeX("$-\\log_{10}(P_{observed})$"),
            key_cols=c("class", "Pvalue"),
            aes_ribbon = aes(ymin=clower, ymax=cupper),
            width=170,
            height=120,
            by_chr=FALSE,
            print_p=FALSE) + facet_wrap(~type)
        print(p)
    }
    dev.off()
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--meta_analysis_results_folder",
    default="~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/meta_results_n_cases_100", required=FALSE,
    help="Folder containing the meta-analysis results files")
parser$add_argument("--out", default="~/Repositories/BRaVa_curation/plots/meta_analysis/meta_analysis_qq_100.pdf", required=FALSE,
    help="Output filepath")
parser$add_argument("--type", default=NULL, required=FALSE,
    help="Which meta-analysis results to plot {'Stouffer', 'weighted Fisher'}")
args <- parser$parse_args()

# Debug
# args <- list()
# args$meta_analysis_results_folder <- "../data/meta_analysis/gcloud"
# args$meta_analysis_results_path_regexp <- "*cleaned*"

main(args)
