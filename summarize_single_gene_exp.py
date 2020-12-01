### Boas Pucker ###
### bp423@cam.ac.uk ###
### v0.2 ###

__usage__ = """
					python summarize_single_gene_exp.py
					--data <COMMA_SEP_LIST_OF_DATA_FILES>
					--species <COMMA_SEP_SPECIES_NAMES>
					--out <OUTPUT_FOLDER>
					
					optional:
					--break <Y_AXIS_BREAKPOINT_PER_GENE>[100]
					--cut <UPPER_CUTOFF_PER_GENE>[10000]
					--genes <COMMA_SEPARATED_GENE_LIST>[CHS,DFR,ANS]
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys, re, glob
import numpy as np
from operator import itemgetter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

# --- end of imports --- #

def load_exp( exp_file ):
	"""! @brief load expression from summary file """
	
	exp = {}
	with open( exp_file, "r" ) as f:
		headers = f.readline().strip().split('\t')[1:]
		for header in headers:
			exp.update( { header: [] } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			for idx, val in enumerate( map( float, parts[1:] ) ):
				exp[ headers[ idx ] ].append( val )
			line = f.readline()
	return exp


def generate_pathway_overview_figure( data, species, figurefile, dataoutputfile, borders, upper_cutoffs, genes ):
	"""! @brief generate boxplots as overview of central flavonoid biosynthesis steps """
		
	fig, axs = plt.subplots( 2, 3, figsize=(12, 5) )
	first_ax_list = axs[0]
	secon_ax_list = axs[1]
	
	gene_labels = [ ]
	for gene in genes:
		gene_labels.append( "$\it{" + gene + "}$" )
	
	with open( dataoutputfile, "w" ) as out:
		out.write( "Gene\t" + "\t".join( species ) + "\n" )
		for idx, gene in enumerate( genes ):
			fax = first_ax_list[ idx ]	#first axis
			sax = secon_ax_list[ idx ]	#second axis
			border = borders[ idx ]	#get the appropriate point where to break per gene
			
			values = []
			positions = []
			upper_cutoff = []
			new_line = []
			
			for k, spec in enumerate( species ):
				positions.append( k )
				values.append( data[ spec ][ gene ] )
				new_line.append( ",".join( map( str, data[ spec ][ gene ] ) ) )
				upper_cutoff.append( 1.1*np.median( data[ spec ][ gene ] ) )
			
			out.write( gene + "\t" + "\t".join( new_line ) + "\n" )
			
			color = "grey"
			sax.boxplot( values, positions=positions, labels=species, widths=0.3,
										patch_artist=True, showfliers=False,	#showmeans=True, notch=True, 
										boxprops=dict(facecolor=color, color=color),
										capprops=dict(color=color),
										whiskerprops=dict(color=color),
										flierprops=dict(color=color, markeredgecolor=color),
										medianprops=dict(color="red")
									)
			fax.boxplot( values, positions=positions, labels=species, widths=0.3,
										patch_artist=True, showfliers=False,	#showmeans=True, notch=True, 
										boxprops=dict(facecolor=color, color=color),
										capprops=dict(color=color),
										whiskerprops=dict(color=color),
										flierprops=dict(color=color, markeredgecolor=color),
										medianprops=dict(color="red")
									)
			
			sax.set_ylabel( gene_labels[ idx ] + " transcript abundance [TPM]", fontsize=8 )
			sax.set_xlim( -0.5, len( species )-0.5 )
			fax.set_xlim( -0.5, len( species )-0.5 )
			
			sax.set_ylim( 0, border )
			fax.set_ylim( border+1, upper_cutoffs[ idx ] )
			
			species_lables = []
			for spec in species:
				species_lables.append( "$\it{" + spec.replace(" ", "}$ $\it{") + "}$" )
			
			sax.set_xticklabels( species_lables, rotation = 90, fontsize=6 )
			
			sax.tick_params(axis='y', which='major', labelsize=8)
			sax.tick_params(axis='y', which='minor', labelsize=8)
			fax.tick_params(axis='y', which='major', labelsize=8)
			fax.tick_params(axis='y', which='minor', labelsize=8)
			
			fax.spines['bottom'].set_visible(False)
		
			sax.spines['top'].set_visible(False)
			fax.set_xticks( [] )
			fax.tick_params(labeltop='off')
			sax.xaxis.tick_bottom()
			
			d = .005  # how big to make the diagonal lines in axes coordinates
			# arguments to pass to plot, just so we don't keep repeating them
			kwargs = dict(transform=fax.transAxes, color='k', clip_on=False)
			fax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
			fax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

			kwargs.update(transform=sax.transAxes)  # switch to the bottom axes
			sax.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
			sax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
		
			
		# legend_elements = []	#species name and sample size
		# for k, spec in enumerate( species ):
			# legend_elements.append( mpatches.Patch( color=colors[ k ], label="$\it{" + spec.replace("_", "}$ $\it{").replace(" ", "}$ $\it{") + "}$" + " (n=" + str( len( all_values[ spec ][0] ) ) + ")" ) )
		
		# ax.legend( handles=legend_elements, loc="upper center", bbox_to_anchor=( 0.5, 1.1 ), ncol=5, fontsize=6 )
		
		plt.subplots_adjust( left=0.05, right=0.99, top=0.99, bottom=0.35, hspace=0.07  )
		
		fig.savefig( figurefile, dpi=600 )
		
		# # --- compare beta vs. all --- #
		# for x, gene in enumerate( genes ):
			# beta_values = data[ "Beta_vulgaris_B2" ][ gene ]
			# for k, spec in enumerate( species ):
				# if k > 0:
					# values = data[ spec ][ gene ]
					# u, p = stats.mannwhitneyu( beta_values, values )
					# print gene + "\tBv vs. " + spec + ": p=" + str( p )


def main( arguments ):
	"""! @brief run everything """
	
	data_files = arguments[ arguments.index('--data')+1 ].split(',')
	species = arguments[ arguments.index('--species')+1 ].split(',')
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if output_folder[-1] != "/":
		output_folder += "/"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	figurefile = output_folder + "TPM_subplots.png"
	dataoutputfile = output_folder + "TPM_subplots.txt"
	
	data = {}	#{ species_name: expression }
	
	
	if len( data_files ) != len( species ):
		sys.exit( "ERROR: number of species names does not match number of provided data files!\nnumber of species names: " + str( len( species ) ) + "\nnumber of data files: " + str( len( data_files ) ) )
	
	for idx, name in enumerate( species ):
		data.update( { name: load_exp( data_files[ idx ] ) } )
	
	if '--break' in arguments:
		borders = map( int, arguments[ arguments.index('--break')+1 ].split(',') )
	else:
		borders = [ 100 ] * len( species )
	
	if '--cut' in arguments:
		upper_cutoffs = map( int, arguments[ arguments.index('--cut')+1 ].split(',') )
	else:
		upper_cutoffs = [ 10000 ] * len( species )
	
	if '--genes' in arguments:
		genes = arguments[ arguments.index('--genes')+1 ].split(',')
	else:
		genes = [ "CHS", "DFR", "ANS" ]
	
	if len( borders ) != len( genes ):
		sys.exit( "ERROR: number of genes does not match number of specified Y axis breaks!\nnumber of genes: " + str( len( genes ) ) + "\nnumber of breaks: " + str( len( borders ) ) )
	
	if len( upper_cutoffs ) != len( genes ):
		sys.exit( "ERROR: number of genes does not match number of specified upper cutoffs!\nnumber of genes: " + str( len( genes ) ) + "\nnumber of upper cutoffs: " + str( len( upper_cutoffs ) ) )
	
	generate_pathway_overview_figure( data, species, figurefile, dataoutputfile, borders, upper_cutoffs, genes )


if '--data' in sys.argv and '--species' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
