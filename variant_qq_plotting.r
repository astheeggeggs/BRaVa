library(data.table)
library(dplyr)
library(ggplot2)

files <- dir("meta_analysis", pattern="rare.exome.tsv.gz", full.names=TRUE)
length <- 10000
maxP <- log10(1)
ribbon_p <- 0.95
source("QC/utils/pretty_plotting.r")

phenotype_class <- list(
	binary = c("AAA","AcApp","AcuLymLeuk","Adenomy","AMD","ALamy","AUD","AloAre",
		"AnoNer","AoSten","Asth","AtopDis","AFib","ADHD","ASD","BCLL","BenCervUterNeo",
		"BenIntNeo","BenNodGoit","BladCanc","BrainCNSCanc","BreastCanc","BrugSynd",
		"BuliNer","BullPemph","CarShock","HCM","CRVO","CervCanc","CML","COPD","CRF",
		"CoffSirSynd","ColonRectCanc","CAD","CCANS","EatDis","Endocar","Endometr",
		"EsophCanc","EssThrom","EFRMB","FSP","FemInf","FemInfAC","FolLymph","Gout",
		"GravesDis","HemoChromo","HF","HepCarcin","HTN","HHD","HypoThyr","HypoThyrSec",
		"IPF","ITP","IBD","IFHern","ILDSarc","IodDef","KabSynd","KidCanc","KleefSynd",
		"LaryxCanc","Leuk","LiverCanc","LiverFibCirr","LongQTSynd","LymphThyrit",
		"MalInf","MatHem","MatHypDis","MODYDiab","MultiMyel","MS","MECS","Myocard",
		"Narco1","NonFuncPitAd","NHL","NonPapTCCBlad","NonRheuValv","OUD","OCD",
		"OvCanc","Pancreat","ParkDis","PeptUlcer","PAD","PlacInsuf","PCOS",
		"PolycythVera","Preeclamps","PregLoss","POAG","PrimSjoSynd","Prolactinom",
		"Psori","RheumHeaDis","RheumArth","RomWardSynd","Sarcoid","SebDerm",
		"SpinaBifAp","StomCanc","Stroke","SLE","TAAD","ThyroCanc","T2Diab","Urolith",
		"UterCanc","VaricVeins","VTE"),
	continuous = c("ALT", "AlcCons", "AST","BMI","CRP","CACS","CK","HDLC","Height",
		"LDLC","TChol","TG","WHRBMI", "LVH","Append","HipRep","CogAbil","EduAtt",
		"PsySymp","SchGrades","SCDCAT")
	)

pdf(file="QQ_variant.pdf", width=6, height=4)
for (file in files)
{
	dt <- fread(file, key="P-value")
	phenotype <- gsub(".*\\/([A-Za-z0-9]+)_.*", "\\1", file)
	
	# We want the QQ to be linear on the log10 scale, so we should sample in that way
	nP <- nrow(dt)
	minP <- log10(1/(nP+1))
	elements <- unique(round(10^seq(minP, maxP, length.out=length) * (nP+1)))
	# Ensure that elements contains the top 100 associations
	elements <- sort(union(elements, seq(1,100)))

	dt_plot <- dt[elements, ]
	dt_plot[ , P_expected := -log10(seq(10^minP, 10^maxP, length.out=(nP+1))[elements])]
	dt_plot <- dt_plot %>% 
		rename(P_observed = `P-value`) %>% 
		mutate(
			P_observed = ifelse(P_observed < 1e-320, 320, -log10(P_observed)),
			clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = (nP:1)[elements], shape1 = (1:nP)[elements])),
			cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = (nP:1)[elements], shape1 = (1:nP)[elements])),
			OR = exp(Effect)
		)

	type <- ifelse(phenotype %in% phenotype_class$continuous, "continuous", "binary")

	if (type == "continuous") {
		dt_plot$color <- cut(dt_plot$Effect,
        breaks = c(-Inf, -0.5, 0, 0.5, Inf),
        labels = c("< -0.5", "[-0.5, 0)", "[0, 0.5]", "> 0.5"))
	    dt_plot$color <- factor(dt_plot$color,
	        levels = c("< -0.5", "[-0.5, 0)", "[0, 0.5]", "> 0.5"))
	    dummy_data <- data.frame(P_expected = NA, P_observed = NA,
	        color = factor(c("< -0.5", "[-0.5, 0)", "[0, 0.5]", "> 0.5"),
	        levels = c("< -0.5", "[-0.5, 0)", "[0, 0.5]", "> 0.5"))
	    )
	} else {
		dt_plot$color <- cut(dt_plot$OR,
        breaks = c(-Inf, 0.5, 1, 2, Inf),
        labels = c("< 0.5", "[0.5, 1)", "[1, 2]", "> 2"))
	    dt_plot$color <- factor(dt_plot$color,
	        levels = c("< 0.5", "[0.5, 1)", "[1, 2]", "> 2"))
	    dummy_data <- data.frame(P_expected = NA, P_observed = NA,
	        color = factor(c("< 0.5", "[0.5, 1)", "[1, 2]", "> 2"),
	        levels = c("< 0.5", "[0.5, 1)", "[1, 2]", "> 2"))
	    )
	}

	p <- create_pretty_qq_plot(
		plot_title=phenotype,
		plot_subtitle="coding regions; gnomAD popmax < 0.01",
		rbind(dt_plot, dummy_data, fill=TRUE),
		aes(x=P_expected, y=P_observed, col=color), key_cols=c("P_observed"),
		aes_ribbon = aes(ymin=clower, ymax=cupper),
		x_label=TeX("$-\\log_{10}(P_{expected})$"), 
		y_label=TeX("$-\\log_{10}(P_{observed})$"),
		print_p=FALSE,
		)

	if (type == "continuous")
	{
		p <- p + scale_color_manual(
            values = c(
                "< -0.5" = "blue3",
                "[-0.5, 0)" = "cornflowerblue",
                "[0, 0.5]" = "indianred3",
                "> 0.5" = "red"),
            labels = c("< -0.5" = "< -0.5",
                "[-0.5, 0)" = "[-0.5, 0)",
                "[0, 0.5]" = "[0, 0.5]",
                "> 0.5" = "> 0.5"),
            name = "Effect size",  aesthetics = c("colour", "fill")
            ) + guides(colour = guide_legend(override.aes = list(size=5)))
	} else {
		p <- p + scale_color_manual(
			values = c(
				"< 0.5" = "blue3",
                "[0.5, 1)" = "cornflowerblue",
                "[1, 2]" = "indianred3",
                "> 2" = "red"),
            labels = c("< 0.5" = "< 0.5",
                "[0.5, 1)" = "[0.5, 1)",
                "[1, 2]" = "[1, 2]",
                "> 2" = "> 2"),
            name = "Odds ratio",  aesthetics = c("colour", "fill")
            ) + guides(colour = guide_legend(override.aes = list(size=5)))
	}
	print(p)
}
dev.off()
