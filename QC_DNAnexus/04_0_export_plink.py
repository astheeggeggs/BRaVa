import hail as hl
import sys
import pyspark
import dxpy
import hail as hl
import pandas as pd
from math import ceil
import json
import os
import re
import glob
import argparse

def export_file(path, out_folder):
    dxpy.upload_local_file(
        filename=path,
        name=path.split('/')[-1],
        folder=out_folder,
        parents=True
    )

def main(chrom):

	my_database = dxpy.find_one_data_object(
		name="my_database", 
		project=dxpy.find_one_project()["id"]
	)["id"]
	database_dir = f'dnax://{my_database}'
	sc = pyspark.SparkContext()
	spark = pyspark.sql.SparkSession(sc)
	hl.init(sc=sc)
     
	def get_final_filter_mt_path(chrom):
			return f'{database_dir}/04_final_filter_v2_write_to_mt/ukb_wes_450k.qced.chr{chrom}.mt'

	# Create sample QC metrics restricted and not restricted (target plus padding) 
	# to the target intervals.

	# Inputs
	MT_HARDCALLS = get_final_filter_mt_path(chrom)
	INITIAL_VARIANT_LIST = f'file:///mnt/project/Barney/qc/02_prefilter_variants/02_variant_qc_chr{chrom}_list.tsv'
	SAMPLE_LIST_INITIAL_QC = f'file:///mnt/project/Barney/qc/03_initial_sample/03_initial_sample_qc_chr{chrom}.tsv'
	HIGH_LD_INTERVALS = sys.argv[4]
	REFERENCE = 'GRCh38'

	# Outputs
	PLINK_FILES = sys.argv[5]

	print("Inputs:")
	print('MT_HARDCALLS; input hard calls matrix table: ', MT_HARDCALLS)
	print('SAMPLE_LIST_INITIAL_QC; target intervals file: ', SAMPLE_LIST_INITIAL_QC)
	print('INITIAL_VARIANT_LIST; padded target intervals file: ', INITIAL_VARIANT_LIST)
	print('HIGH_LD_INTERVALS; set of high-LD regions for removal: ', HIGH_LD_INTERVALS)

	print("Outputs:")
	print('PLINK_FILES; prefix for high-quality common variant plink files: ', PLINK_FILES)

	hl.init(default_reference=REFERENCE)

	ht_initial_samples = hl.import_table(SAMPLE_LIST_INITIAL_QC, no_header=True, key='f0')
	ht_initial_variants = hl.import_table(INITIAL_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome=REFERENCE), 'alleles':hl.tarray(hl.tstr)})
	ht_initial_variants = ht_initial_variants.key_by(ht_initial_variants.locus, ht_initial_variants.alleles)

	high_LD_intervals = hl.import_locus_intervals(HIGH_LD_INTERVALS, reference_genome=REFERENCE)

	mt = hl.read_matrix_table(MT_HARDCALLS)
	mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
	mt = mt.annotate_rows(in_high_LD = hl.is_defined(high_LD_intervals[mt.locus]))

	mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]) & (~mt.in_high_LD))
	mt = mt.filter_rows(mt.locus.in_x_nonpar() | mt.locus.in_autosome_or_par())
	mt = hl.variant_qc(mt, name='qc')
	mt = mt.filter_rows(
		(mt.qc.AF[0] > 0.01) & 
		(mt.qc.AF [0]< 0.99) & 
		((mt.qc.call_rate > 0.98) | mt.locus.in_x_nonpar() | mt.locus.in_x_par())
		).persist()

	mt.count()

	mt_chr = hl.filter_intervals(
		mt, [hl.parse_locus_interval(hl.eval('chr' + hl.str(x)), reference_genome=REFERENCE)])
	n_chr = mt_chr.count_rows()

	print('\nn variants in chr')
	print(x)
	print(n_chr)

	hl.export_plink(mt_chr, PLINK_FILES + '.chr' + str(x))
