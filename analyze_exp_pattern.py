### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.35 ###

__usage__ = """
					python analyze_exp_pattern.py
					--in <SUMMARY_FILE>
					--out <OUTPUT_FOLDER>
					--taxon <TAXON_FILE>
					--genes <COMMA_SEPARATED_GENE_LIST>
					"""

import re, os, sys, math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
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


def generate_pathway_overview_figure_broken_yaxis( O, A, B, figurefile, O_group, A_group, B_group, genes, border, max_y_cutoff=2000 ):
	"""! @brief generate boxplots as overview of entire flavonoid biosynthesis pathway """
	
	universal_fontsize = 14
	sample_size_fontsize = 8
	
	fig, ( ax, ax2 ) = plt.subplots( 2, 1, figsize=(10, 5) )	#, sharex=True
	
	O_positions = []
	A_positions = []
	B_positions = []
	O_values = []
	A_values = []
	B_values = []
	labels = []
	
	print figurefile
	for idx, gene in enumerate( genes ):
		labels.append( "$\it{" + gene.replace( "-", "'") +"}$"  )
		O_positions.append( idx )
		A_positions.append( idx+0.25 )
		B_positions.append( idx+0.5 )
		O_values.append( O[ gene ] )
		A_values.append( A[ gene ] )
		B_values.append( B[ gene ] )
		
		try:
			print gene
			
			try:
				print "O mean : " + str( np.median( O[ gene ] ) )
			except:
				pass
			try:
				print "A mean : " + str( np.median( A[ gene ] ) )
			except:
				pass
			try:
				print "B mean : " + str( np.median( B[ gene ] ) )
			except:
				pass
			
			r, p = stats.mannwhitneyu( O[ gene ], A[ gene ] )
			print "outgroup vs. group A: " + str( p )
			r, p = stats.mannwhitneyu( O[ gene ], B[ gene ] )
			print "outgroup vs. group B: " + str( p )
			r, p = stats.mannwhitneyu( A[ gene ], B[ gene ] )
			print "group A vs. group B: " + str( p )
			print ""
		except ValueError:
			pass
	
	O_color = "grey"
	ax.boxplot( O_values, positions=O_positions, widths=0.2,
						patch_artist=True, showfliers=False,	#showmeans=True, notch=True, 
						boxprops=dict(facecolor=O_color, color=O_color),
						capprops=dict(color=O_color),
						whiskerprops=dict(color=O_color),
						flierprops=dict(color=O_color, markeredgecolor=O_color),
						medianprops=dict(color="black"),
						meanprops=dict(color="black"),
					)
	ax2.boxplot( O_values, positions=O_positions, widths=0.2,
						patch_artist=True, showfliers=False,	#showmeans=True, notch=True, 
						boxprops=dict(facecolor=O_color, color=O_color),
						capprops=dict(color=O_color),
						whiskerprops=dict(color=O_color),
						flierprops=dict(color=O_color, markeredgecolor=O_color),
						medianprops=dict(color="black"),
						meanprops=dict(color="black"),
					)
	B_color = "magenta"
	ax.boxplot( B_values, positions=B_positions, widths=0.2,
						patch_artist=True, showfliers=False,	#showmeans=True, notch=True, 
						boxprops=dict(facecolor=B_color, color=B_color),
						capprops=dict(color=B_color),
						whiskerprops=dict(color=B_color),
						flierprops=dict(color=B_color, markeredgecolor=B_color),
						medianprops=dict(color="black") 
					)
	ax2.boxplot( B_values, positions=B_positions, widths=0.2,
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
	ax2.boxplot( A_values, positions=A_positions, labels=labels, widths=0.2,				#plotting this block last places labels in center
						patch_artist=True, showfliers=False,	#showmeans=True, notch=True, 
						boxprops=dict(facecolor=A_color, color=A_color),
						capprops=dict(color=A_color),
						whiskerprops=dict(color=A_color),
						flierprops=dict(color=A_color, markeredgecolor=A_color),
						medianprops=dict(color="black") 
					)
	
	# # --- adding sample sizes --- #
	# all_Ys = []
	# for idx, sample in enumerate( A_values ):
		# n = len( sample )
		# if n > 0:
			# y = np.percentile( sample, 75 )
		# else:
			# y = 0
		# x = A_positions[ idx ]
		# if n > 0:
			# if y >= border:
				# ax.text( x, y + 0.01*max_y_cutoff, "n="+str( n ), fontsize=sample_size_fontsize, rotation=90, va="bottom" )
			# else:
				# ax2.text( x, y+ 0.01*border, "n="+str( n ), fontsize=sample_size_fontsize, rotation=90, va="bottom" )
		# all_Ys.append( y )
	
	# for idx, sample in enumerate( B_values ):
		# n = len( sample )
		# if n > 0:
			# y = np.percentile( sample, 75 )
		# else:
			# y = 0
		# x = B_positions[ idx ]
		# if n > 0:
			# if y >= border:
				# ax.text( x, y + 0.01*max_y_cutoff, "n="+str( n ), fontsize=sample_size_fontsize, rotation=90, va="bottom" )
			# else:
				# ax2.text( x, y+ 0.01*border, "n="+str( n ), fontsize=sample_size_fontsize, rotation=90, va="bottom" )
		# all_Ys.append( y )
	
	# for idx, sample in enumerate( O_values ):
		# n = len( sample )
		# if n > 0:
			# y = np.percentile( sample, 75 )
		# else:
			# y = 0
		# x = O_positions[ idx ]
		# if n > 0:
			# if y >= border:
				# ax.text( x, y + 0.01*max_y_cutoff, "n="+str( n ), fontsize=sample_size_fontsize, rotation=90, va="bottom" )
			# else:
				# ax2.text( x, y+ 0.01*border, "n="+str( n ), fontsize=sample_size_fontsize, rotation=90, va="bottom" )
		# all_Ys.append( y )
	
	ax.set_ylabel( "gene expression [TPMs]", fontsize=universal_fontsize )
	ax.set_xlim( -0.3, len( genes )-0.2 )
	ax2.set_xlim( -0.3, len( genes )-0.2 )
	
	ax.tick_params(axis='both', which='major', labelsize=universal_fontsize )
	ax.tick_params(axis='both', which='minor', labelsize=universal_fontsize )
	ax2.tick_params(axis='both', which='major', labelsize=universal_fontsize )
	ax2.tick_params(axis='both', which='minor', labelsize=universal_fontsize )
	
	
	
	ax2.set_ylim( 0, border )
	ax.set_ylim( border+1, max_y_cutoff )	#max( [ x for sublist in A_values + B_values + O_values for x in sublist ] )
	
	
	ax.spines['bottom'].set_visible(False)
	
	ax2.spines['top'].set_visible(False)
	ax.set_xticks( [] )
	ax.tick_params(labeltop='off')
	ax2.xaxis.tick_bottom()
	
	d = .005  # how big to make the diagonal lines in axes coordinates
	# arguments to pass to plot, just so we don't keep repeating them
	kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
	ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
	ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

	kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
	ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
	ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
	
	if len( [ x for sublist in O_values for x in sublist ] ) > 0:
		legend_elements = [ 	
											mpatches.Patch( color=O_color, label= "outgroup"  ),
											mpatches.Patch( color=A_color, label= "anthocyanin" ),	#A_group
											mpatches.Patch( color=B_color, label= "betalain" ) 	#B_group
										]
	else:
		legend_elements = [ 	
											mpatches.Patch( color=A_color, label= "anthocyanin" ),	#A_group
											mpatches.Patch( color=B_color, label= "betalain" ) 	#B_group
										]
	ax.legend( handles=legend_elements, loc="upper center", bbox_to_anchor=( 0.5, 0.995 ), ncol=3, fontsize=universal_fontsize )
	
	plt.subplots_adjust( left=0.1, right=0.999, top=0.98, bottom=0.06, hspace=0.03 )
	
	fig.savefig( figurefile )


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
	generate_pathway_overview_figure_broken_yaxis( normed_data_O, normed_data_A, normed_data_B, output_folder + "pathway_overview_A2_B2.png", "A0", "A2", "B2", genes, 100 )
	
	normed_data_O, normed_data_A, normed_data_B = generate_data_matrix_file( data, taxon_mapping, output_folder + "table_for_heatmap_A0_A2_B3.txt", "A0", "A2", "B3", genes )
	generate_pathway_overview_figure_broken_yaxis( normed_data_O, normed_data_A, normed_data_B, output_folder + "pathway_overview_A2_B3.png", "A0", "A2", "B3", genes, 100  )
	
	
	normed_data_O, normed_data_A, normed_data_B = generate_data_matrix_file( data, taxon_mapping, output_folder + "table_for_heatmap_A0_A2_B4.txt", "A0", "A2", "B4", genes )
	generate_pathway_overview_figure_broken_yaxis( normed_data_O, normed_data_A, normed_data_B, output_folder + "pathway_overview_A2_B4.png", "A0", "A2", "B4", genes, 100  )
	
	
	normed_data_O, normed_data_A, normed_data_B = generate_data_matrix_file( data, taxon_mapping, output_folder + "table_for_heatmap_X_A_B.txt", "Y", "A", "B", genes )
	generate_pathway_overview_figure_broken_yaxis( normed_data_O, normed_data_A, normed_data_B, output_folder + "pathway_overview_X_A_B.png", "Y", "A", "B", genes, 100  )
	


if '--in' in sys.argv and '--out' in sys.argv and '--taxon' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
