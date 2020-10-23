### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python coexp_analysis.py
					--exp <EXPRESSION_FILE>
					--refgene <REF_GENE_ID(s)>
					--gene1 <GENE1_ID(s)>
					--gene2 <GENE2_ID(s)>
					--gene3 <GENE3_ID(s)>
					--out <OUTPUT_FILE>
					
					optional:
					--refname <REF_GENE_DISPLAY_NAME>
					--name1 <GENE1_DISPLAY_NAME>
					--name2 <GENE2_DISPLAY_NAME>
					--name3 <GENE3_DISPLAY_NAME>
					"""

import math, sys, os
import scipy.stats as stats

# --- end of imports --- #


def load_exp( exp_file ):
	"""! @brief load expression from summary file """
	
	exp = {}
	with open( exp_file, "r" ) as f:
		first_row = f.readline().strip()
		if "gene" in first_row:
			headers = first_row.split('\t')[1:]
		else:
			headers = first_row.split('\t')
		for header in headers:
			exp.update( { header: {} } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			for idx, val in enumerate( map( float, parts[1:] ) ):
				exp[ headers[ idx ] ].update( { parts[0]: val+0.001 } )
			line = f.readline()
	return exp


def generate_output_file( ref_values1, ref_values2, ref_values3, gene1_values, gene2_values, gene3_values, ref_name, name1, name2, name3, output_file ):
	"""! @brief generate figure """
	
	if len( ref_values1 ) > 0:
		rs_gene1, ps_gene1 = stats.spearmanr( ref_values1, gene1_values )
		rp_gene1, pp_gene1 = stats.pearsonr( ref_values1, gene1_values )
	else:
		rs_gene1, ps_gene1 = "nan", "nan"
		rp_gene1, pp_gene1 = "nan", "nan"
	
	if len( ref_values2 ) > 0:
		rs_gene2, ps_gene2 = stats.spearmanr( ref_values2, gene2_values )
		rp_gene2, pp_gene2 = stats.pearsonr( ref_values2, gene2_values )
	else:
		rs_gene2, ps_gene2 = "nan", "nan"
		rp_gene2, pp_gene2 = "nan", "nan"
	
	if len( ref_values3 ) > 0:
		rs_gene3, ps_gene3 = stats.spearmanr( ref_values3, gene3_values )
		rp_gene3, pp_gene3 = stats.pearsonr( ref_values3, gene3_values )
	else:
		rs_gene3, ps_gene3 = "nan", "nan"
		rp_gene3, pp_gene3 = "nan", "nan"
	
	
	with open( output_file, "w" ) as out:
		out.write( "\t".join( [ "RefGene", "Gene", "r_Spearman", "p_Spearman", "r_Pearson", "p_Pearson", "SampleSize" ] ) + "\n" )
		out.write( "\t".join( map( str, [ ref_name, name1, rs_gene1, ps_gene1, rp_gene1, pp_gene1, len( ref_values1 ) ] ) ) + "\n" )
		out.write( "\t".join( map( str, [ ref_name, name2, rs_gene2, ps_gene2, rp_gene2, pp_gene2, len( ref_values2 ) ] ) ) + "\n" )
		out.write( "\t".join( map( str, [ ref_name, name3, rs_gene3, ps_gene3, rp_gene3, pp_gene3, len( ref_values3 ) ] ) ) + "\n" )


def main( arguments ):
	"""! @brief run everything """
	
	tpm_file = arguments[ arguments.index( '--exp' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	
	if '/' in arguments[ arguments.index( '--refgene' )+1 ]:
		ref_genes = arguments[ arguments.index( '--refgene' )+1 ].split('/')
	else:
		ref_genes = [ arguments[ arguments.index( '--refgene' )+1 ] ]
	
	if '/' in arguments[ arguments.index( '--gene1' )+1 ]:
		gene1 = arguments[ arguments.index( '--gene1' )+1 ].split('/')
	else:
		gene1 = [ arguments[ arguments.index( '--gene1' )+1 ] ]
	if '/' in arguments[ arguments.index( '--gene2' )+1 ]:
		gene2 = arguments[ arguments.index( '--gene2' )+1 ].split('/')
	else:
		gene2 = [ arguments[ arguments.index( '--gene2' )+1 ] ]
	if '/' in arguments[ arguments.index( '--gene3' )+1 ]:
		gene3 = arguments[ arguments.index( '--gene3' )+1 ].split('/')
	else:
		gene3 = [ arguments[ arguments.index( '--gene3' )+1 ] ]
	
	
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
	if '--name3' in arguments:
		name3 = arguments[ arguments.index( '--name3' )+1 ]
	else:
		name3 = "Gene3"
	
	
	exp_data = load_exp( tpm_file )
	ref_values1 = []
	ref_values2 = []
	ref_values3 = []
	gene1_values = []
	gene2_values = []
	gene3_values = []

	for sample in exp_data.keys():
		ref_gene_counter = 0
		for gene in ref_genes:
			try:
				ref_gene_counter += exp_data[ sample ][ gene ]
			except KeyError:
				pass	#print sample + "\t" + gene
		gene1_counter = 0
		for gene in gene1:
			try:
				gene1_counter += exp_data[ sample ][ gene ]
			except KeyError:
				pass	#print sample + "\t" + gene
		
		gene2_counter = 0
		for gene in gene2:
			try:
				gene2_counter += exp_data[ sample ][ gene ]
			except KeyError:
				pass	#print sample + "\t" + gene
		
		gene3_counter = 0
		for gene in gene3:
			try:
				gene3_counter += exp_data[ sample ][ gene ]
			except KeyError:
				pass	#print sample + "\t" + gene
		
		if ref_gene_counter > 0 and gene1_counter > 0:
			ref_values1.append( math.log( ref_gene_counter ) )
			gene1_values.append( math.log( gene1_counter ) )
		
		if ref_gene_counter > 0 and gene2_counter > 0:
			ref_values2.append( math.log( ref_gene_counter ) )
			gene2_values.append( math.log( gene2_counter ) )
		
		if ref_gene_counter > 0 and gene3_counter > 0:
			ref_values3.append( math.log( ref_gene_counter ) )
			gene3_values.append( math.log( gene3_counter ) )

	generate_output_file( ref_values1, ref_values2, ref_values3, gene1_values, gene2_values, gene3_values, ref_name, name1, name2, name3, output_file )


if '--exp' in sys.argv and '--out' in sys.argv and '--refgene' in sys.argv and '--gene1' in sys.argv and '--gene2' in sys.argv and '--gene3' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
