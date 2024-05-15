#!/usr/bin/env bash
#SBATCH --job-name=meta-analysis-qq
#SBATCH --output=/well/lindgren/dpalmer/logs/%x_%j.log
#SBATCH --error=/well/lindgren/dpalmer/logs/%x_%j.err

META_FILE_GENE=$1
OUT=$2

Rscript meta_analysis_qq.r --meta_analysis_results_folder $META_FILE_GENE --out $OUT --include_gene_names
