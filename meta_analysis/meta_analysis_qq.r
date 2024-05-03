#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(latex2exp)
library(argparse)

source("~/Repositories/BRaVa_curation/QC/utils/pretty_plotting.r")
gene_name_mapping_file <- "~/Repositories/BRaVa_curation/data/gene_mapping.txt.gz"
T <- 6 # -log10(P) threshold for inclusion on the plots

main <- function(args)
{
    ribbon_p <- 0.95
    files <- dir(args$meta_analysis_results_folder, full.names=TRUE)
    files_to_include <- grep("meta", dir(args$meta_analysis_results_folder))
    files <- files[files_to_include]
    print(files)
    print(args)

    if (args$include_gene_names) {
        if (file.exists(gene_name_mapping_file)) {
            dt_genes <- fread(gene_name_mapping_file)
            if (!(all(c("Gene stable ID", "Gene name") %in% names(dt_genes)))) {
                warning(paste0("gene name mapping file does not contain 'Gene stable ID' ",
                    "(ensembl ID *without* version) and 'Gene name' (gene symbol for plotting), "
                    "reverting to not plitting gene-names"))
                 args$include_gene_names <- FALSE
            } else {
                names(dt_genes)[which(names(dt_genes) == "Gene stable ID")] <- "Region"
                names(dt_genes)[which(names(dt_genes) == "Gene name")] <- "labels"
                setkey(dt_genes, "Region")
            }
        } else {
            warning("gene name mapping file is not present, reverting to not plotting gene-names")
            args$include_gene_names <- FALSE
        }
    }

    pdf(file=args$out, width=8, height=8)
    for (file in files) {
        phe <- gsub(".*/(.*)_meta.*", "\\1", file)
        phe_plot <- gsub("_", " ", gsub("_$", "", str_trim(gsub("[[:space:]_]+", "\\_", phe))))
        cat(paste0(phe_plot, "...\n"))
        dt_meta <- fread(cmd = paste("gzcat", file), key="Region")
        dt_meta[, Pvalue := -log10(Pvalue)]
        if(args$include_gene_names) {
            dt_meta <- merge(dt_meta, dt_genes, all.x=TRUE)
        }
        
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
                variant_class_plot <- gsub("_", " ", gsub("[\\|;]", ",\n", g))
                cex_labels <- 2
                dt_current <- dt_meta_to_plot %>% filter(Group==g, max_MAF==m)
                p <- create_pretty_qq_plot(
                    plot_title=phe_plot,
                    plot_subtitle=paste0(variant_class_plot, "; max MAF = ", TeX(max_MAF_plot)),
                    cex_labels=cex_labels,
                    dt_current,
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
                
                if(args$include_gene_names) {
                    p <- p + geom_label_repel(data=dt_current %>% filter(class == "Burden", Pvalue > T),
                        aes(label=labels), box.padding = 0.5, label.padding=0.1, point.padding = 0.2,
                        color = 'grey30', segment.color = 'grey50',
                        size=cex_labels, segment.size=0.1, show.legend = FALSE)
                }   
                p <- p + facet_wrap(~type)
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
            print_p=FALSE)

        if(args$include_gene_names) {
            p <- p + geom_label_repel(data=dt %>% filter(class == "Burden", Pvalue > T),
                aes(label=labels), box.padding = 0.5, label.padding=0.1, point.padding = 0.2,
                color = 'grey30', segment.color = 'grey50',
                size=cex_labels, segment.size=0.1, show.legend = FALSE)
        } 
        p <- p + facet_wrap(~type)
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
parser$add_argument("--include_gene_names", default=FALSE, action='store_true',
    help="Do you want to highlight the most significant genes with their gene-names?")
args <- parser$parse_args()

# Debug
# args <- list()
# args$meta_analysis_results_folder <- "../data/meta_analysis/gcloud"
# args$meta_analysis_results_path_regexp <- "*cleaned*"

main(args)
