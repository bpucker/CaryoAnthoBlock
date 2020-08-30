### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python check_gene_exp_correlation.py
					--exp <EXPRESSION_FILE>
					--refgene <REF_GENE_ID(s)>
					--gene1 <GENE1_ID(s)>
					--gene2 <GENE2_ID(s)>
					--out <OUTPUT_FILE>
					
					optional:
					--refname <REF_GENE_DISPLAY_NAME>
					--name1 <GENE1_DISPLAY_NAME>
					--name2 <GENE2_DISPLAY_NAME>
					"""


import matplotlib.pyplot as plt
import math, sys, os
import scipy.stats as stats
import matplotlib.patches as mpatches

# --- end of imports --- #


def load_exp( exp_file ):
	"""! @brief load expression from summary file """
	
	exp = {}
	with open( exp_file, "r" ) as f:
		headers = f.readline().strip().split('\t')[1:]
		for header in headers:
			exp.update( { header: {} } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			for idx, val in enumerate( map( float, parts[1:] ) ):
				if val > 0:
					exp[ headers[ idx ] ].update( { parts[0]: math.log( val ) } )
				else:
					exp[ headers[ idx ] ].update( { parts[0]: val } )
			line = f.readline()
	return exp


def generate_figure( ref_values, gene1_values, gene2_values, ref_name, name1, name2, output_file ):
	"""! @brief generate figure """
	
	fig, ax = plt.subplots()

	ax.scatter( ref_values, gene1_values, s=10, color="orange" )
	ax.scatter( ref_values, gene2_values, s=10, color="dodgerblue" )

	r_gene1, p_gene1 = stats.spearmanr( ref_values, gene1_values )
	r_gene2, p_gene2 = stats.spearmanr( ref_values, gene2_values )
	
	ax.text( min( ref_values ), 0.9*( max( gene1_values+gene2_values ) ), ref_name+" vs. "+name1+": r=" + str(  r_gene1)[:4] + " (p=" + str( round( p_gene1, 5 ) ) + ")", fontsize=10 )
	ax.text( min( ref_values ), 0.8*( max( gene1_values+gene2_values ) ), ref_name+" vs. "+name2+": r=" + str(  r_gene2)[:4] + " (p=" + str( round( p_gene2, 5 ) ) + ")", fontsize=10 )

	ax.set_xlabel( ref_name + " transcript abundance [log() of TPMs]" )
	ax.set_ylabel( name1 + "/"+ name2 +" transcript abundance [log() of TPMs]" )
	
	my_legend = [ mpatches.Patch(color="orange", label=name1 ), mpatches.Patch(color="dodgerblue", label=name2 ) ]
	ax.legend( handles=my_legend, loc="upper center", ncol=2, bbox_to_anchor=( 0.5, 1.1 ), fontsize=10 )

	fig.savefig( output_file, dpi=300 )


def main( arguments ):
	"""! @brief run everything """
	
	tpm_file = arguments[ arguments.index( '--exp' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	if ',' in arguments[ arguments.index( '--refgene' )+1 ]:
		ref_gene = arguments[ arguments.index( '--refgene' )+1 ].split(',')
	else:
		ref_gene = [ arguments[ arguments.index( '--refgene' )+1 ] ]
	if ',' in arguments[ arguments.index( '--gene1' )+1 ]:
		gene1 = arguments[ arguments.index( '--gene1' )+1 ].split(',')
	else:
		gene1 = [ arguments[ arguments.index( '--gene1' )+1 ] ]
	if ',' in arguments[ arguments.index( '--gene2' )+1 ]:
		gene2 = arguments[ arguments.index( '--gene2' )+1 ].split(',')
	else:
		gene2 = [ arguments[ arguments.index( '--gene2' )+1 ] ]
	if '--refname' in arguments:
		ref_name = arguments[ arguments.index( '--refname' )+1 ]
	else:
		ref_name = "RefGene"
	if '--name1' in arguments:
		name1 = arguments[ arguments.index( '--name1' )+1 ]
	else:
		name1 = "Gene1"
	if '--name2' in arguments:
		name2 = arguments[ arguments.index( '--name2' )+1 ]
	else:
		name2 = "Gene2"
	
	
	exp_data = load_exp( tpm_file )
	ref_values = []
	gene1_values = []
	gene2_values = []

	for sample in exp_data.keys():
		ref_gene_counter = 0
		for gene in ref_gene:
			try:
				ref_gene_counter += exp_data[ sample ][ gene ]
			except KeyError:
				print sampe + "\t" + gene
		gene1_counter = 0
		for gene in gene1:
			try:
				gene1_counter += exp_data[ sample ][ gene ]
			except KeyError:
				print sampe + "\t" + gene
		
		gene2_counter = 0
		for gene in gene2:
			try:
				gene2_counter += exp_data[ sample ][ gene ]
			except KeyError:
				print sampe + "\t" + gene
		
		ref_values.append( ref_gene_counter )
		gene1_values.append( gene1_counter )
		gene2_values.append( gene2_counter )

	generate_figure( ref_values, gene1_values, gene2_values, ref_name, name1, name2, output_file )


if '--exp' in sys.argv and '--out' in sys.argv and '--refgene' in sys.argv and '--gene1' in sys.argv and '--gene2' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
