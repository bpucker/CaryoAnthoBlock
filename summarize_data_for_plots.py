### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.33 ###

__usage__ = """
					python summarize_data_for_plots.py
					--in <SUMMARY_FILE>
					--out <OUTPUT_FOLDER>
					--taxon <TAXON_FILE>
					--genes <COMMA_SEPARATED_GENE_LIST>
					
					optional:
					--ymax <Y_AXIS_CUTOFF>
					"""

import re, os, sys, math
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import matplotlib.patches as mpatches

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


def generate_data_matrix_file( data, taxon_mapping, output_file, A_group, B_group, genes ):
	"""! @brief generate a data matrix for heatmap construction """
	
	normed_data_A = {}	#A lineage
	normed_data_B = {}	#B lineage
	with open( output_file, "w" ) as out:
		for gene in genes:
			normed_data_A.update( { gene: [] } )
			normed_data_B.update( { gene: [] } )
		out.write( "\t".join( [ "spec" ] + genes ) + "\n" )
		for species in data.keys():
			try:
				spec = taxon_mapping[ species ]
			except KeyError:
				spec = "XXX"
			
			if spec[:len( A_group )] == A_group or spec[:len( B_group )] == B_group:
				values_for_new_line = []
				for gene in genes:
					try:
						values_for_new_line.append( data[ species ][ gene ]['median'] )
					except KeyError:
						values_for_new_line.append( 0 )
				new_line = [ spec ]
				for idx, val in enumerate( values_for_new_line ):
					new_line.append( val )
					if spec[:len( A_group )] == A_group:
						if val > -1:
							normed_data_A[ genes[ idx ] ].append( val )
					elif spec[:len( B_group )] == B_group:
						if val > -1:
							normed_data_B[ genes[ idx ] ].append( val )
					else:
						print "ERROR: value was not assigned to any group: " + species + "\t" + genes[ idx ]
				out.write( "\t".join( map( str, new_line ) ) + "\n" )
				
	return normed_data_A, normed_data_B


def generate_figure_and_table( A, B, figure_file, dataoutputfile, A_group, B_group, genes, cutoff ):
	"""! @brief generate boxplots as overview of entire flavonoid biosynthesis pathway """
	
	# --- generation of data output file --- #
	with open( dataoutputfile, "w" ) as out:
		out.write( "Gene\tAnthocyanins\tBetalains\n" )
		for idx, gene in enumerate( genes ):
			new_line = [ gene, ",".join( map( str, A[ gene ] ) ), ",".join( map( str, B[ gene ] ) ) ]
			out.write( "\t".join( new_line ) + "\n" )
	
	# --- generation of figures --- #
	universal_fontsize = 14
	sample_size_fontsize = 8
	
	fig, ax = plt.subplots( figsize=(10, 5) )	#, sharex=True
	
	A_positions = []
	B_positions = []
	A_values = []
	B_values = []
	labels = []
	
	for idx, gene in enumerate( genes ):
		labels.append( "$\it{" + gene.replace( "-", "'") +"}$"  )
		A_positions.append( idx )
		B_positions.append( idx+0.25 )
		A_values.append( A[ gene ] )
		B_values.append( B[ gene ] )
		
		try:
			print gene
			try:
				print "A mean : " + str( np.median( A[ gene ] ) )
			except:
				pass
			try:
				print "B mean : " + str( np.median( B[ gene ] ) )
			except:
				pass
			
			r, p = stats.mannwhitneyu( A[ gene ], B[ gene ] )
			print A_group + " vs. " +  B_group + ": " + str( p ) + "\n"
		except ValueError:
			pass
	
	B_color = "magenta"
	ax.boxplot( B_values, positions=B_positions, widths=0.2,
						patch_artist=True, showfliers=False,	#showmeans=True, notch=True, 
						boxprops=dict(facecolor=B_color, color=B_color),
						capprops=dict(color=B_color),
						whiskerprops=dict(color=B_color),
						flierprops=dict(color=B_color, markeredgecolor=B_color),
						medianprops=dict(color="black") 
					)
	
	A_color = "dodgerblue"
	ax.boxplot( A_values, positions=A_positions, labels=labels, widths=0.2,				#plotting this block last places labels in center
						patch_artist=True, showfliers=False,	#showmeans=True, notch=True, 
						boxprops=dict(facecolor=A_color, color=A_color),
						capprops=dict(color=A_color),
						whiskerprops=dict(color=A_color),
						flierprops=dict(color=A_color, markeredgecolor=A_color),
						medianprops=dict(color="black") 
					)

	ax.set_ylabel( "gene expression [TPMs]", fontsize=universal_fontsize )
	ax.set_xlim( -0.3, len( genes )-0.2 )
	ax.set_ylim( 0, cutoff )
	
	ax.tick_params(axis='both', which='major', labelsize=universal_fontsize )
	ax.tick_params(axis='both', which='minor', labelsize=universal_fontsize )
	
	legend_elements = [ 	mpatches.Patch( color=A_color, label= "anthocyanin" ),	#A_group
										mpatches.Patch( color=B_color, label= "betalain" ) 	#B_group
									]
	ax.legend( handles=legend_elements, loc="upper right", bbox_to_anchor=( 0.999, 0.999 ), ncol=1, fontsize=universal_fontsize )
	
	plt.subplots_adjust( left=0.1, right=0.999, top=0.98, bottom=0.06, hspace=0.03 )
	
	fig.savefig( figure_file )


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index( '--in' )+1 ]
	output_folder = arguments[ arguments.index( '--out' )+1 ]
	taxon_table = arguments[ arguments.index( '--taxon' )+1 ]
	genes = arguments[ arguments.index( '--genes' )+1 ].split(',')
	
	if '--ymax' in arguments:
		cutoff = int( arguments[ arguments.index( '--ymax' )+1 ] )
	else:
		cutoff = 10000
	
	data = load_data_from_table( input_file )
	origin = True
	taxon_mapping = load_spec_name_mapping_table( taxon_table, origin )
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- generate data table for heatmap construction --- #
	normed_data_A, normed_data_B = generate_data_matrix_file( data, taxon_mapping, output_folder + "table_for_heatmap_A_B2.txt", "A", "B2", genes )
	generate_figure_and_table( normed_data_A, normed_data_B, output_folder + "pw_A_B2.pdf", output_folder + "pw_A_B2.data.txt", "A", "B2", genes, cutoff )
	
	normed_data_A, normed_data_B = generate_data_matrix_file( data, taxon_mapping, output_folder + "table_for_heatmap_A_B3.txt", "A", "B3", genes )
	generate_figure_and_table( normed_data_A, normed_data_B, output_folder + "pw_A_B3.pdf", output_folder + "pw_A_B3.data.txt", "A", "B3", genes, cutoff  )
	
	
	normed_data_A, normed_data_B = generate_data_matrix_file( data, taxon_mapping, output_folder + "table_for_heatmap_A_B4.txt", "A", "B4", genes )
	generate_figure_and_table( normed_data_A, normed_data_B, output_folder + "pw_A_B4.pdf", output_folder + "pw_A_B4.data.txt", "A", "B4", genes, cutoff  )
	
	
	normed_data_A, normed_data_B = generate_data_matrix_file( data, taxon_mapping, output_folder + "table_for_heatmap_A_B.txt", "A", "B", genes )
	generate_figure_and_table( normed_data_A, normed_data_B, output_folder + "pw_A_B.pdf", output_folder + "pw_A_B.data.txt", "A", "B", genes, cutoff  )
	


if '--in' in sys.argv and '--out' in sys.argv and '--taxon' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
