### Boas Pucker ###
### bp423@cam.ac.uk ###
### v0.2 ###

__usage__ = """
					python prep_TPM_summary_pairwise.py
					--exp1 <BV_TPM_FILE>
					--config1 <BV_CONFIG_FILE>
					--exp2 <AT_TPM_FILE>
					--config2 <ATH_CONFIG_FILE>
					--out <OUTPUT_FOLDER>
					
					optional:
					--spec1 <SPECIES_NAME1; _ will be replaced by space>
					--spec2 <SPECIES_NAME2; _ will be replaced by space>
					--genes <COMMA_SEPARATED_LIST_OF_GENE_IDs>[flavonoid biosynthesis genes]
										
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys, re, glob, math
import numpy as np
from operator import itemgetter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

# --- end of imports --- #


def load_config( config_file ):
	"""! @brief load pathway assignment and order of genes, filenames with homologs in all species """
	
	gene_order = []
	candidate_genes = {}
	with open( config_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				try:
					genes = parts[2]
					if "," in genes:
						genes = genes.split(',')
					else:
						genes = [ genes ]
					for gene in genes:
						candidate_genes.update( { gene: parts[1]} )
				except:
					print line
			line = f.readline()
	return gene_order, candidate_genes


def load_bv_expression_per_seq( counttable, candidate_genes ):
	"""! @brief load TPMs per sequence of given species count table """
	
	exp = {}
	with open( counttable, "r" ) as f:
		samples = f.readline().strip().split('\t')	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				name = candidate_genes[ parts[0] ]
				for idx, val in enumerate( parts[1:] ):
					if not "nan" in val:
						val = float( val )
						try:
							exp[ name ][ samples[ idx ] ] += val
						except KeyError:
							try:
								exp[ name ].update( { samples[ idx ]: math.log( val+0.0001 ) } )
							except KeyError:
								exp.update( { name: { samples[ idx ]: math.log( val+0.0001 ) } } )
			except KeyError:
				pass
			
			line = f.readline()
	return exp, samples


def load_ath_expression_per_seq( counttable, candidate_genes ):
	"""! @brief load TPMs per sequence of given species count table """
	
	exp = {}
	with open( counttable, "r" ) as f:
		samples = f.readline().strip().split('\t')[1:]	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				name = candidate_genes[ parts[0] ]
				for idx, val in enumerate( parts[1:] ):
					if not "nan" in val:
						val = float( val )
						try:
							exp[ name ][ samples[ idx ] ] += val
						except KeyError:
							try:
								exp[ name ].update( { samples[ idx ]: math.log( val+0.0001 ) } )
							except KeyError:
								exp.update( { name: { samples[ idx ]: math.log( val+0.0001 ) } } )
			except KeyError:
				pass
			line = f.readline()
	return exp, samples


def generate_data_matrix_file( data, samples, output_file, genes ):
	"""! @brief generate a data matrix for heatmap construction """
	
	normed_data = {}	#B lineage
	with open( output_file, "w" ) as out:
		for gene in genes:
			normed_data.update( { gene: [] } )
		out.write( "\t".join( [ "spec" ] + genes ) + "\n" )
		for sample in samples:
			values_for_new_line = []
			for gene in genes:
				try:
					values_for_new_line.append( data[ gene ][ sample ] )
				except KeyError:
					values_for_new_line.append( 0 )
			cum_tpm = float( sum( values_for_new_line ) )
			new_line = [ sample ]
			for idx, val in enumerate( values_for_new_line ):
				new_line.append( val )
				normed_data[ genes[ idx ] ].append( val )
			out.write( "\t".join( map( str, new_line ) ) + "\n" )
	return normed_data


def generate_pathway_overview_figure( bv_normed_data, at_normed_data, figurefile, species1, species2, genes, capping_max_value=800, border=40 ):
	"""! @brief generate boxplots as overview of entire flavonoid biosynthesis pathway """
	
	fig, ax = plt.subplots( figsize=(10, 5) )
	
	bv_positions = []
	bv_values = []
	at_positions = []
	at_values = []
	labels = []
	
	for idx, gene in enumerate( genes ):
		labels.append( "$\it{" + gene.replace( "-", "'") +"}$"  )
		at_positions.append( idx )
		bv_positions.append( idx+0.5 )
		bv_values.append( bv_normed_data[ gene ] )
		at_values.append( at_normed_data[ gene ] )
	
		try:
			print gene
			r, p = stats.mannwhitneyu( bv_normed_data[ gene ], at_normed_data[ gene ] )
			print species1 + " vs. " + species2 + ": " + str( p )
			print np.median( at_normed_data[ gene ] ) / np.median( bv_normed_data[ gene ] )
			print ""
		except:
			pass
	
	at_color = "dodgerblue"
	bv_color = "magenta"
	
	ax.boxplot( at_values, positions=at_positions, labels=labels, widths=0.3,
						notch=True, patch_artist=True, showfliers=False,	#showmeans=True, 
						boxprops=dict(facecolor=at_color, color=at_color),
						capprops=dict(color=at_color),
						whiskerprops=dict(color=at_color),
						flierprops=dict(color=at_color, markeredgecolor=at_color),
						medianprops=dict(color="black") 
					)
	ax.boxplot( bv_values, positions=bv_positions, labels=labels, widths=0.3,				#plotting this block last places labels in center
						notch=True, patch_artist=True, showfliers=False,	#showmeans=True, 
						boxprops=dict(facecolor=bv_color, color=bv_color),
						capprops=dict(color=bv_color),
						whiskerprops=dict(color=bv_color),
						flierprops=dict(color=bv_color, markeredgecolor=bv_color),
						medianprops=dict(color="black") 
					)
	
	ax.set_ylabel( "gene expression [log10() of TPM]", fontsize=16 )
	ax.tick_params(axis='both', which='major', labelsize=16)
	ax.tick_params(axis='both', which='minor', labelsize=16)
	
	ax.set_xlim( -0.3, len( genes )-0.2 )
	ax.set_ylim( -10, 50 )

	legend_elements = [ 	mpatches.Patch( color=at_color, label="$\it{" + species2.replace("_", "}$ $\it{") + "}$" + " (n=" + str( len( at_values[0] ) ) + ")" ),
									mpatches.Patch( color=bv_color, label="$\it{" + species1.replace("_", "}$ $\it{") + "}$" + " (n=" + str( len( bv_values[0] ) ) + ")" )								
									]
	
	ax.legend( handles=legend_elements, loc="upper center", bbox_to_anchor=( 0.5, 0.995 ), ncol=2, fontsize=16 )
	
	plt.subplots_adjust( left=0.15, right=0.99, top=0.97, bottom=0.075, hspace=0.03  )
	fig.savefig( figurefile )


def main( arguments ):
	"""! @brief run everything """
	
	bv_exp_file = arguments[ arguments.index('--exp1')+1 ]
	bv_config_file = arguments[ arguments.index('--config1')+1 ]
	at_exp_file = arguments[ arguments.index('--exp2')+1 ]
	at_config_file = arguments[ arguments.index('--config2')+1 ]
	
	if '--spec1' in arguments:
		species1 = arguments[ arguments.index('--spec1')+1 ]
	else:
		species1 = "species1"
	if '--spec2' in arguments:
		species2 = arguments[ arguments.index('--spec2')+1 ]
	else:
		species2 = "species2"
	
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if '--genes' in arguments:
		genes = arguments[ arguments.index('--out')+1 ]
		if "," in genes:
			genes = genes.split(',')
		else:
			genes = [ genes ]
	else:
		genes = [ "PAL", "C4H", "4CL", "CHS", "CHI", "F3H", "F3-H", "F3-5-H", "DFR", "ANS", "LAR", "ANR" ]
	
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- loading all information from config file --- #
	bv_gene_order, bv_candidate_genes = load_config( bv_config_file )
	at_gene_order, at_candidate_genes = load_config( at_config_file )
	
	# --- load expression data --- #
	bv_exp_data, bv_samples = load_bv_expression_per_seq( bv_exp_file, bv_candidate_genes )
	at_exp_data, at_samples = load_ath_expression_per_seq( at_exp_file, at_candidate_genes )
	
	# --- generate data matrix file --- #
	bv_output_file = output_folder + species1 + "_data_matrix_for_heatmap.txt"
	bv_normed_data = generate_data_matrix_file( bv_exp_data, bv_samples, bv_output_file, genes )
	
	at_output_file = output_folder + species2 + "_data_matrix_for_heatmap.txt"
	at_normed_data = generate_data_matrix_file( at_exp_data, at_samples, at_output_file, genes ) 
	
	# --- generate boxplot overview figure --- #
	overview_boxplot_figure = output_folder + species1.replace( " ", "_" ) + "_vs_" + species2.replace( " ", "_" ) + "_boxplot.png"
	generate_pathway_overview_figure( bv_normed_data, at_normed_data, overview_boxplot_figure, species1, species2, genes )


if '--exp1' in sys.argv and '--config1' in sys.argv and '--out' in sys.argv and '--config2' in sys.argv and '--exp2' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
