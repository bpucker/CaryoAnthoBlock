### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.31 ###

__usage__ = """
					python summarize_data_for_plots.py
					--in <SUMMARY_FILE>
					--out <OUTPUT_FOLDER>
					--taxon <TAXON_FILE>
					--genes <COMMA_SEPARATED_GENE_LIST>
					"""

import re, os, sys, math
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# --- end of imports --- #

def load_data_from_table( input_file ):
	""""! @brief load data from given table """
	
	data = {}
	with open( input_file, "r" ) as f:
		headers = f.readline().strip().split('\t')[2:]
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			data.update( { parts[0]: {} } )
			for idx, part in enumerate( parts[2:] ):
				gene = headers[ idx ].split('_')[0]
				if idx % 4 == 0:	#median
					data[ parts[0] ].update( { gene: { 'median': float( part ) } } )
				elif idx % 4 == 1:	#IQR
					data[ parts[0] ][ gene ].update( { 'IQR': float( part ) } )
				elif idx % 4 == 2:	#copy number
					data[ parts[0] ][ gene ].update( { 'copies': float( part ) } )
				elif idx % 4 == 3:	#single values
					if len( part.strip() ) > 0:
						if ";" in part:
							data[ parts[0] ][ gene ].update( { 'values': map( float, part.split(';') ) } )
						else:
							data[ parts[0] ][ gene ].update( { 'values': [ float( part ) ]  } )
					else:
						pass	#print part
			line = f.readline()
	return data


def load_spec_name_mapping_table( taxon_file, origin ):
	"""! @brief load species name mapping table (include origin in name if specified) """
	
	mapping_table = {}
	with open( taxon_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if origin:
				mapping_table.update( { parts[1]: parts[3] + "_" + parts[2] } )
			else:
				mapping_table.update( { parts[1]: parts[2] } )
			line = f.readline()
	return mapping_table


def generate_data_matrix_file( data, taxon_mapping, output_file, O_group, A_group, B_group, genes ):
	"""! @brief generate a data matrix for heatmap construction """
	
	normed_data_O = {}	#outgroup
	normed_data_A = {}	#A lineage
	normed_data_B = {}	#B lineage
	with open( output_file, "w" ) as out:
		
		for gene in genes:
			normed_data_O.update( { gene: [] } )
			normed_data_A.update( { gene: [] } )
			normed_data_B.update( { gene: [] } )
		out.write( "\t".join( [ "spec" ] + genes ) + "\n" )
		for species in data.keys():
			try:
				spec = taxon_mapping[ species ]
			except KeyError:
				spec = "XXX"
			
			if spec[:len(O_group  )] in [ O_group, A_group, B_group]:
				values_for_new_line = []
				for gene in genes:
					try:
						values_for_new_line.append( data[ species ][ gene ]['median'] )
					except KeyError:
						values_for_new_line.append( 0 )
				new_line = [ spec ]
				for idx, val in enumerate( values_for_new_line ):
					new_line.append( val )
					if spec[:len(O_group  )] == O_group:
						if val > -1:
							normed_data_O[ genes[ idx ] ].append( val )
					elif spec[:len(O_group  )] == A_group:
						if val > -1:
							normed_data_A[ genes[ idx ] ].append( val )
					else:
						if val > -1:
							normed_data_B[ genes[ idx ] ].append( val )
				out.write( "\t".join( map( str, new_line ) ) + "\n" )
				
	return normed_data_O,normed_data_A, normed_data_B


def generate_summary_table( O, A, B, dataoutputfile, O_group, A_group, B_group, genes ):
	"""! @brief generate boxplots as overview of entire flavonoid biosynthesis pathway """
	
	with open( dataoutputfile, "w" ) as out:
		out.write( "Gene\tOutgroup\tAnthocyanins\tBetalains\n" )
		for idx, gene in enumerate( genes ):
			new_line = [ gene, ",".join( map( str, O[ gene ] ) ), ",".join( map( str, A[ gene ] ) ), ",".join( map( str, B[ gene ] ) ) ]
			out.write( "\t".join( new_line ) + "\n" )


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index( '--in' )+1 ]
	output_folder = arguments[ arguments.index( '--out' )+1 ]
	taxon_table = arguments[ arguments.index( '--taxon' )+1 ]
	genes = arguments[ arguments.index( '--genes' )+1 ].split(',')
	
	data = load_data_from_table( input_file )
	origin = True
	taxon_mapping = load_spec_name_mapping_table( taxon_table, origin )
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- generate data table for heatmap construction --- #
	normed_data_O, normed_data_A, normed_data_B = generate_data_matrix_file( data, taxon_mapping, output_folder + "table_for_heatmap_A0_A2_B2.txt", "A0", "A2", "B2", genes )
	generate_summary_table( normed_data_O, normed_data_A, normed_data_B, output_folder + "pathway_overview_A2_B2.data.txt", "A0", "A2", "B2", genes )
	
	normed_data_O, normed_data_A, normed_data_B = generate_data_matrix_file( data, taxon_mapping, output_folder + "table_for_heatmap_A0_A2_B3.txt", "A0", "A2", "B3", genes )
	generate_summary_table( normed_data_O, normed_data_A, normed_data_B, output_folder + "pathway_overview_A2_B3.data.txt", "A0", "A2", "B3", genes  )
	
	
	normed_data_O, normed_data_A, normed_data_B = generate_data_matrix_file( data, taxon_mapping, output_folder + "table_for_heatmap_A0_A2_B4.txt", "A0", "A2", "B4", genes )
	generate_summary_table( normed_data_O, normed_data_A, normed_data_B, output_folder + "pathway_overview_A2_B4.data.txt", "A0", "A2", "B4", genes  )
	
	
	normed_data_O, normed_data_A, normed_data_B = generate_data_matrix_file( data, taxon_mapping, output_folder + "table_for_heatmap_X_A_B.txt", "Y", "A", "B", genes )
	generate_summary_table( normed_data_O, normed_data_A, normed_data_B, output_folder + "pathway_overview_X_A_B.data.txt", "Y", "A", "B", genes  )
	


if '--in' in sys.argv and '--out' in sys.argv and '--taxon' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
