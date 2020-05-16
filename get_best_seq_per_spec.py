### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python get_best_seq_per_spec.py
					--in <INPUT_FOLDER>
					--out <OUTPUT_FOLDER>
					
					optional:
					--ratio <FLOAT, MINIMAL_PROPORTION_OF_CONSERVED_RESIDUES>[0]
					--cutoff <INT, MINIMAL_NUMBER_OF_CONSERVED_RESIDUES>[-1]
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys, glob
from operator import itemgetter

# --- end of imports --- #

def get_best_sequence_per_species( filename, ratio, cutoff ):
	"""! @brief produce list with best matching sequence per species """
	
	best_hit = {}
	best_matches = []	#list of best matching sequences (max one per species if missmatches)
	with open( filename, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 1:
				matches = parts[1:].count( "True" )
				if matches == len( parts[1:] ):
					best_matches.append( parts[0] )
				spec = parts[0].split('@')[0]
				try:
					if matches > best_hit[ spec ]['matches']:
						best_hit[ spec ] = { 'ID': parts[0], 'matches': matches }
				except KeyError:
					if matches > cutoff and matches / float( len( parts[1:] ) ) > ratio:
						best_hit.update( { spec: { 'ID': parts[0], 'matches': matches } } )
			line = f.readline()
	for each in best_hit.values():
		best_matches.append( each['ID'] )
	return list( set( best_matches ) )


def main( arguments ):
	"""! @brief run everything """
	
	input_folder = arguments[ arguments.index('--in')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if '--ratio' in arguments:
		ratio = float( arguments[ arguments.index('--ratio')+1 ] )
	else:
		ratio = 0
	
	if '--cutoff' in arguments:
		cutoff = int( arguments[ arguments.index('--cutoff')+1 ] )
	else:
		cutoff = 0
	
	if output_folder[-1] != '/':
		output_folder += "/"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	input_files = glob.glob( input_folder + "*.txt" )
	for filename in input_files:
		ID = filename.split('/')[-1].split('.')[0]
		best_seq_per_spec = get_best_sequence_per_species( filename, ratio, cutoff )
		output_file = output_folder + ID + ".txt"
		with open( output_file, "w" ) as out:
			out.write( "\n".join( best_seq_per_spec )  )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
