import hail as hl
import sys

# python variant_extraction_for_plotting.py ${RESULTS_FILE} ${RESULTS_FILE_RARE_EXOME}

# Inputs
REFERENCE = 'GRCh38'
PROTEIN_CODING_REGIONS = "data/protein_coding_genes.txt.gz"
GNOMAD_POPMAX = "data/gnomad.exomes.r2.1.1.sites.liftover_grch38_popmax_0.01.tsv.bgz"
RESULTS_FILE = sys.argv[1] #'meta_analysis/FemInf_variant_meta_analysis1.tbl'

# Outputs
PROTEIN_CODING_INTERVALS = 'data/protein_coding_genes_intervals.txt.bgz'
RESULTS_FILE_RARE_EXOME = sys.argv[2] #'meta_analysis/FemInf_variant_meta_analysis1.rare.exome.tsv.gz'

hl.init(default_reference=REFERENCE)
# export JAVA_HOME="$(/usr/libexec/java_home -v 1.8)"

ht = hl.import_table(PROTEIN_CODING_REGIONS, impute=True, delimiter=' ', force=True)

# Define the list of vald chromosomes
valid_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']

# Filter the table to keep only the rows with valid chromosomes
ht = ht.filter(hl.literal(valid_chromosomes).contains(ht.chromosome_name))

# Create a set of coding intervals to filter the variant results to
ht = ht.annotate(contig = 'chr' + ht.chromosome_name)
ht = ht.select(ht.contig, ht.start_position, ht.end_position)
ht.export(PROTEIN_CODING_INTERVALS, header=False)

ht_results = hl.import_table(RESULTS_FILE, impute=True, delimiter='\t')
ht_results = ht_results.filter(ht_results.Allele1 != 'ur')
ht_results = ht_results.annotate(variant = hl.parse_variant(ht_results['MarkerName']))

coding_intervals = hl.import_locus_intervals(PROTEIN_CODING_INTERVALS)
ht_results = ht_results.annotate(in_coding_interval = hl.is_defined(coding_intervals[ht_results.variant.locus]))
ht_results = ht_results.filter(ht_results.in_coding_interval)

# Finally, remove all variants that exceed 1% in gnomAD
ht_gnomAD_popmax = hl.import_table(GNOMAD_POPMAX, no_header=True)
ht_gnomAD_popmax = ht_gnomAD_popmax.annotate(variant = hl.parse_variant(ht_gnomAD_popmax.f0))
ht_gnomAD_popmax = ht_gnomAD_popmax.key_by('variant').select()

ht_results = ht_results.annotate(in_gnomAD_popmax = hl.is_defined(ht_gnomAD_popmax[ht_results.variant]))
ht_results = ht_results.filter(~ht_results.in_gnomAD_popmax)
ht_results = ht_results.key_by('P-value')
ht_results = ht_results.select(
	ht_results.MarkerName,
	ht_results.Allele1,
	ht_results.Allele2, ht_results.Effect,
	ht_results.StdErr, ht_results.Direction,
	ht_results.HetISq, ht_results.HetChiSq, 
	ht_results.HetDf, ht_results.HetPVal)
ht_results.export(RESULTS_FILE_RARE_EXOME)
