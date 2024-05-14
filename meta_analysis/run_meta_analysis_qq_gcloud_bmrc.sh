#!/usr/bin/env bash
#$ -N meta-analysis
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/
#$ -e /well/lindgren/dpalmer/logs/
#$ -P lindgren.prjc
#$ -q short.qe@@short.hge

FILES_GENE=$1
OUT=$2

Rscript meta_analysis.r --file_paths $FILES_GENE --no_sex_check --out $OUT
