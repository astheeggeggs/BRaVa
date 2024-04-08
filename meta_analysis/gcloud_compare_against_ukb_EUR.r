library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(latex2exp)

source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")

data_directory <- "~/Repositories/BRaVa_curation/data/meta_analysis/gcloud"

chr <- c(seq(1,22), "X")
pops <- "EUR"
max_MAF_groups <- c(1e-4, 1e-3, 1e-2)
type_plot <- "weighted Fisher"
groups <- c(
    "damaging_missense_or_protein_altering", 
    "other_missense_or_protein_altering",
    "pLoF",
    "pLoF;damaging_missense_or_protein_altering",
    "pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
    "synonymous"
)
case_control_threshold <- 100
tests <- c("Burden", "SKAT", "SKAT-O")
phenotype <- "Coronary_artery_disease"
file <- "uk-biobank.palmer.PRELIMINARY.Coronary_artery_disease.JULY23Freeze.ALL.EUR.19915.382460.SAIGE.gene.20240110.txt.gz"

pdf(file=paste0("~/Repositories/BRaVa_curation/plots/meta_analysis/gcloud_meta_analysis_EUR_comparison_", case_control_threshold, "_cutoff.pdf"),
    width=6, height=6)
dt <- fread(paste0(data_directory, "/", file)) %>% mutate(max_MAF = as.character(max_MAF))
dt_annots <- dt %>% filter(Group != "Cauchy")
dt_cauchy <- dt %>% filter(Group == "Cauchy")
p_trim <- gsub("_$", "", str_trim(gsub("[[:space:]_]+", "\\_", phenotype)))
dt_meta <- fread(cmd=paste0(
    "gzcat ", data_directory,
    paste0(
        "/meta_results/n_cases_", case_control_threshold, "/",
        p_trim, "_gene_meta_analysis_", case_control_threshold,
        "_cutoff.tsv.gz")))
if (!is.null(type_plot)) {
    dt_meta <- dt_meta %>% filter(type == type_plot)
}
dt_annots <- melt(
    dt_annots, measure.vars=c("Pvalue_Burden", "Pvalue_SKAT", "Pvalue"),
        variable.name = "class", value.name = "Pvalue_EUR")
dt_annots <- dt_annots %>% mutate(
    class = ifelse(class == "Pvalue", "SKAT-O",
            ifelse(class == "Pvalue_SKAT", "SKAT", "Burden")
        )
    )
setkeyv(dt_annots, c("Region", "Group", "max_MAF", "class"))
setkeyv(dt_meta, c("Region", "Group", "max_MAF", "class"))
dt_merge <- merge(dt_annots, dt_meta)

for (m in max_MAF_groups) {
    for (g in groups) {
        pl <- ggplot(dt_merge %>% filter(Group == g, max_MAF == m),
            aes(x=-log10(Pvalue_EUR), y=-log10(Pvalue))) + 
        geom_bin_2d(bins=100) + geom_abline(intercept=0,slope=1,col='red',lwd=0.3) +
        facet_grid(rows=vars(class), cols=vars(type)) +
        theme_bw() +
        labs(
            x=TeX("$-\\log_{10}(P_{EUR})$"),
            y=TeX("$-\\log_{10}(P)$"),
            title=gsub("_", " ", p_trim),
            subtitle=g
            )
        print(pl)
    }
}

dt_cauchy <- melt(dt_cauchy, measure.vars=c("Pvalue_Burden", "Pvalue_SKAT", "Pvalue"),
        variable.name = "class", value.name = "Pvalue_EUR")
dt_cauchy <- dt_cauchy %>% mutate(
    class = ifelse(class == "Pvalue", "SKAT-O",
        ifelse(class == "Pvalue_SKAT", "SKAT", "Burden")
        )) %>% filter(Group == "Cauchy") %>% select(-c("max_MAF", "Group"))
dt_meta <- dt_meta %>% filter(Group == "Cauchy") %>% select(-c("max_MAF", "Group"))
setkeyv(dt_cauchy, c("Region", "class"))
setkeyv(dt_meta, c("Region", "class"))
dt_merge <- merge(dt_cauchy, dt_meta)

pl <- ggplot(dt_merge,
    aes(x=-log10(Pvalue_EUR), y=-log10(Pvalue))) + 
geom_bin_2d(bins=100) + geom_abline(intercept=0,slope=1,col='red',lwd=0.3) +
facet_grid(rows=vars(class), cols=vars(type)) +
theme_bw() +
labs(
    x=TeX("$-\\log_{10}(P_{EUR})$"),
    y=TeX("$-\\log_{10}(P)$"),
    title=gsub("_", " ", p_trim),
    subtitle="Cauchy"
    )
print(pl)
dev.off()
