### Boas Pucker ###
### bp423@cam.ac.uk ###
### v0.2 ###

__usage__ = """
					python TPM_assess_whole_pathway.py
					--exp <FOLDER_WITH_TPM_FILES>
					--config <CONFIG_FILE>
					--out <OUTPUT_FOLDER>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys, re, glob
import numpy as np
from operator import itemgetter

# --- end of imports --- #

def load_gene_IDs( filename ):
	"""! @brief load all gene IDs from given file """
	
	with open( filename, "r" ) as f:
		content = f.read()
	genes = re.findall( "[a-zA-Z0-9]+@\d+", content )
	return genes


def load_config( config_file ):
	"""! @brief load pathway assignment and order of genes, filenames with homologs in all species """
	
	config = []
	candidate_genes = []
	with open( config_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				genes = load_gene_IDs( parts[2] )
				config.append( { 'pathway': parts[0], 'name': parts[1], 'genes': genes, 'order': len( config ) } )
				candidate_genes.append( parts[1] )
			line = f.readline()
	return config, candidate_genes


def load_expression_per_seq( counttable ):
	"""! @brief load TPMs per sequence of given species count table """
	
	exp = {}
	IQR = {}
	with open( counttable, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			values = map( float, parts[1:] )
			exp.update( { parts[0]: np.median( values ) } )
			IQR.update( { parts[0]: np.percentile( values, 75, interpolation = 'midpoint' ) - np.percentile( values, 25, interpolation = 'midpoint' ) } )
			line = f.readline()
	return exp, IQR, len( parts[1:] )


def get_TPMs_of_seq_of_interest( tpms, IQR, config ):
	"""! @brief get expression and IQR per sequence of interest """
	
	exp_results = {}
	IQR_results = {}
	copy_number_results = {}
	all_exp_results = {}
	for entry in config:
		tmp_exp = []
		tmp_IQR = []
		for gene in entry['genes']:
			try:
				tmp_exp.append( tpms[ gene ] )
				tmp_IQR.append( IQR[ gene ] )
			except KeyError:
				pass
		if len( tmp_exp ) == 1:
			exp_results.update( { entry['name']: tmp_exp[0] } )
			IQR_results.update( { entry['name']: tmp_IQR[0] } )
			copy_number_results.update( { entry['name']: 1 } )
		elif len( tmp_exp ) > 1:
			exp_results.update( { entry['name']: sum( tmp_exp ) } )
			IQR_results.update( { entry['name']: np.median( tmp_IQR ) } )	#alternative: use largest IQR
			copy_number_results.update( { entry['name']: len( tmp_exp ) } )
		else:
			exp_results.update( { entry['name']: -1 } )
			IQR_results.update( { entry['name']: -1 } )
			copy_number_results.update( { entry['name']: -1 } )
		all_exp_results.update( { entry['name']: tmp_exp } )
	return exp_results, IQR_results, copy_number_results, all_exp_results


def generate_summary_output_file( species_list, exp_data_of_interest, exp_data_of_interest_IQRs, copy_number_per_gene_per_spec, single_exp_values_per_spec, summary_table_file, candidate_genes, sample_sizes_per_spec ):
	"""! @brief generate summary output file """
	
	genes = []
	for gene in candidate_genes:	#each gene name appears twice (TPM + IQR)
		genes.append( gene + "_median" )
		genes.append( gene + "_IQR" )
		genes.append( gene + "_CopyNumber" )
		genes.append( gene + "_SingleExpValues" )
	with open( summary_table_file, "w" ) as out:
		out.write( "Species\tNumberOfSamples\t" + "\t".join( genes ) + '\n' )
		for species in sorted( species_list ):
			new_line = [ species, sample_sizes_per_spec[species] ]
			for gene in candidate_genes:
				new_line.append( exp_data_of_interest[ species ][ gene ] )
				new_line.append( exp_data_of_interest_IQRs[ species ][ gene ] )
				new_line.append( copy_number_per_gene_per_spec[ species ][ gene ] )
				if len( single_exp_values_per_spec[ species ][ gene ] ) > 0:
					new_line.append( ";".join( map( str, single_exp_values_per_spec[ species ][ gene ] ) ) )
				else:
					new_line.append( "-1" )
			out.write( "\t".join( map( str, new_line ) ) + '\n' )


def find_dual_pigment_candidates( species_list, exp_data_of_interest, exp_data_of_interest_IQRs, copy_number_per_gene_per_spec, dual_pigment_candidate_spec_file, candidate_genes, sample_sizes_per_spec ):
	"""! @brief find dual pigment candidates """
	
	# --- prepare values for ranking of gene expression values --- #
	values_per_gene = {}
	for gene in candidate_genes:
		values_per_gene.update( { gene: [] } )
	for spec in exp_data_of_interest.keys():
		for gene in candidate_genes:
			try:
				values_per_gene[ gene ].append( exp_data_of_interest[ spec ][ gene ]  )
			except KeyError:
				pass
	
	sorted_values_per_gene = {}
	for gene in candidate_genes:
		sorted_values_per_gene.update( { gene: sorted( values_per_gene[ gene ] )[::-1] } )	#values sorted decreasing; small index = high value
	
	# --- calculate rank for each gene in each species --- #
	ranked_gene_exp_per_spec = {}
	for species in species_list:
		tmp_ranks = {}
		for gene in candidate_genes:
			index = sorted_values_per_gene[ gene ].index( exp_data_of_interest[ species ][ gene ] )
			tmp_ranks.update( { gene: index } )
		ranked_gene_exp_per_spec.update( { species: tmp_ranks } )
	
	# --- sort species by max or average rank --- #
	data_for_sorting = []
	for spec in ranked_gene_exp_per_spec.keys():
		ranks = ranked_gene_exp_per_spec[ spec ].values()
		data_for_sorting.append( { 'spec': spec, 'max': max( ranks ), 'mean': np.mean( ranks ), 'median': np.median( ranks ) } )
	sorted_spec = sorted( data_for_sorting, key=itemgetter('max') )
	
	genes = []
	for gene in candidate_genes:
		genes.append( gene + "_median" )
		genes.append( gene + "_IQR" )
	
	with open( dual_pigment_candidate_spec_file, "w" ) as out:
		out.write( "\t".join( [ "Species", "NumerOfSamples", "MaxRank", "MeanRank", "MedianRank" ] + genes ) + '\n' )
		for spec in sorted_spec:
			new_line = [ spec['spec'], sample_sizes_per_spec[ spec['spec'] ], spec['max'], spec['mean'], spec['median'] ]
			for gene in candidate_genes:
				new_line.append( exp_data_of_interest[ spec['spec'] ][ gene ] )
				new_line.append( exp_data_of_interest_IQRs[ spec['spec']][ gene ] )
			out.write( "\t".join( map( str, new_line ) ) + '\n' )


def main( arguments ):
	"""! @brief run everything """
	
	exp_file_folder = arguments[ arguments.index('--exp')+1 ]
	config_file = arguments[ arguments.index('--config')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- loading all information from config file --- #
	config, candidate_genes = load_config( config_file )
	print "number of candidate genes: " + str( len( candidate_genes ) )
	count_tables_per_transcriptome = glob.glob( exp_file_folder + "*.txt" )
	print "number of detected count tables: " + str( len( count_tables_per_transcriptome ) )
	
	# --- prepare dictionaries to store data of all genes and all species --- #
	exp_data_of_interest = {}
	exp_data_of_interest_IQRs = {}
	copy_number_per_gene_per_spec = {}
	sample_sizes_per_spec = {}
	single_exp_values_per_spec = {}
	#one key per gene
	#cumulative value of all sequences per species => calculated in get_TPMs_of_seq_of_interest( tpms, IQR, config )
	#average value over samples => calculated in load_expression_per_seq( counttable )
	
	# --- load data from individual count tables --- #
	species_list = []
	for k, counttable in enumerate( count_tables_per_transcriptome ):
		print "progress: " + str( k+1 ) + "/" + str( len( count_tables_per_transcriptome ) )
		ID = counttable.split('/')[-1].split('.')[0]
		species_list.append( ID )
		tpms, IQR, sample_size = load_expression_per_seq( counttable )
		sample_sizes_per_spec.update( { ID: sample_size } )
		seq_of_interest_tpms, seq_of_interest_IQRs, copy_number_results, all_exp_results = get_TPMs_of_seq_of_interest( tpms, IQR, config )
		for gene in candidate_genes:
			try:
				exp_data_of_interest[ ID ]
				try:
					exp_data_of_interest[ ID ].update( { gene: seq_of_interest_tpms[ gene ] } )
					exp_data_of_interest_IQRs[ID].update( { gene: seq_of_interest_IQRs[ gene ] } )
					copy_number_per_gene_per_spec[ ID ].update( { gene: copy_number_results[ gene ] } )
					single_exp_values_per_spec[ ID ].update( { gene: all_exp_results[ gene ] } )
				except KeyError:
					exp_data_of_interest[ID].update( { gene: 0 } )
					exp_data_of_interest_IQRs[ID].update( { gene: 0 } )
					copy_number_per_gene_per_spec[ ID ].update( { gene: 0 } )
					single_exp_values_per_spec[ ID ].update( { gene: [] } )
			except KeyError:
				try:
					exp_data_of_interest.update( { ID: { gene: seq_of_interest_tpms[ gene ] } } )
					exp_data_of_interest_IQRs.update( { ID: { gene: seq_of_interest_IQRs[ gene ] } } )
					copy_number_per_gene_per_spec.update( { ID : { gene: copy_number_results[ gene ] } } )
					single_exp_values_per_spec.update( { ID : { gene: all_exp_results[ gene ] } } )
				except KeyError:
					exp_data_of_interest.update( { ID: { gene: 0 } } )
					exp_data_of_interest_IQRs.update( { ID: { gene: 0 } } )
					copy_number_per_gene_per_spec.update( { ID : { gene: 0 } } )
					single_exp_values_per_spec.update( { ID : { gene: [] } } )
	
	# --- generate summary output file --- #
	summary_table_file = output_folder + "summary_table.txt"
	generate_summary_output_file( species_list, exp_data_of_interest, exp_data_of_interest_IQRs, copy_number_per_gene_per_spec, single_exp_values_per_spec, summary_table_file, candidate_genes, sample_sizes_per_spec )
	
	# --- identify dual pigment candidate species --- #
	dual_pigment_candidate_spec_file = output_folder + "dual_pigment_candidates.txt"
	find_dual_pigment_candidates( species_list, exp_data_of_interest, exp_data_of_interest_IQRs, copy_number_per_gene_per_spec, dual_pigment_candidate_spec_file, candidate_genes, sample_sizes_per_spec )


if '--exp' in sys.argv and '--config' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
