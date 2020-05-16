### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python sample_classifier.py
					--pep <TARGET_PEPTIDE_FILE>
					--exp <TARGET_EXPRESSION_FILE>
					--ref <REFERENCE_PEPTIDE_FILE>
					--out <OUTPUT_FOLDER>
					
					optional:
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys, re, glob
import numpy as np
from operator import itemgetter

# --- end of imports --- #


def load_best_blast_hit( blast_result_file ):
	"""! @brief load best blast hit per query """
	
	best_hits = {}
	
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				data = best_hits[ parts[0] ]
				if float( parts[-1] ) > data['score']:
					del best_hits[ parts[0] ]
					best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'ID': parts[1] } } )
			except:
				best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'ID': parts[1] } } )
			line = f.readline()
	
	final = {}
	for key in best_hits.keys():
		final.update( { key: best_hits[ key ]['ID'] } )
	return final


def group_sequence_IDs_by_tissue( best_hits ):
	"""! @brief group sequence IDs by tissue """
	
	groups = {}
	for key in best_hits.keys():
		ID = key.split('@')[0]
		try:
			groups[ ID ].append( best_hits[ key ] )
		except KeyError:
			groups.update( { ID: [ best_hits[ key ] ] } )
	return groups


def load_expression_per_sample( target_exp_file ):
	"""! @brief load expression values per sample """
	
	data = {}
	with open( target_exp_file, "r" ) as f:
		headers = f.readline().strip().split('\t')
		for header in headers:
			data.update( { header: {} } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			values = map( float, parts[1:] )
			for idx, value in enumerate( values ):
				data[ headers[ idx ] ].update( { parts[0]: value } )
			line = f.readline()
	return data


def identify_tissue_per_sample( seq_IDs_per_tissue, exp ):
	"""! @brief identify most likely tissue type per sample """
	
	classification = {}
	for sample in exp.keys():
		exp_per_gene = exp[ sample ]
		candidates = []
		for key in seq_IDs_per_tissue.keys():
			seq_IDs = seq_IDs_per_tissue[ key ]
			tmp = []
			for seq in seq_IDs:
				tmp.append( exp_per_gene[ seq ] )
			candidates.append( { 'tissue': key, 'mean': np.mean( tmp ), 'median': np.median( tmp ), 'max': max( tmp ) } )
		tissue = sorted( candidates, key=itemgetter('mean') )[-1]['tissue']
		classification.update( { sample: tissue } )
	return classification


def main( arguments ):
	"""! @brief run everything """
	
	target_pep_file = arguments[ arguments.index('--pep')+1 ]
	target_exp_file = arguments[ arguments.index('--exp')+1 ]
	reference_pep_file = arguments[ arguments.index('--ref')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if output_folder[-1] != '/':
		output_folder += "/"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	
	blastdb = output_folder + "blastdb"
	blast_result_file = output_folder + "BLAST_results.txt"
	if not os.path.isfile( blast_result_file ):
		os.popen( "makeblastdb -in " + target_pep_file + " -out " + blastdb + " -dbtype prot" )
		os.popen( "blastp -query " + reference_pep_file + " -db " + blastdb + " -out " + blast_result_file + " -outfmt 6 -evalue 0.0001" )
	best_hits = load_best_blast_hit( blast_result_file )
	seq_IDs_per_tissue = group_sequence_IDs_by_tissue( best_hits )
	
	exp = load_expression_per_sample( target_exp_file )
	
	tissue_per_sample = identify_tissue_per_sample( seq_IDs_per_tissue, exp )
	print tissue_per_sample
	output_file = output_folder + "SUMMARY_RESULTS.txt"
	with open( output_file, "w" ) as out:
		out.write( "ID\tTissue\n" )
		for key in tissue_per_sample.keys():
			out.write( key + "\t" + tissue_per_sample[ key ] + '\n' )




if '--pep' in sys.argv and '--exp' in sys.argv and '--ref' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
