### Boas Pucker ###
### bp423@cam.ac.uk ###
### v0.1 ###

__usage__ = """
					python merge_counttables_per_assembly.py
					--in <INPUT_FOLDER>
					--tpms <TPM_OUTPUT_FOLDER>
					--counts <COUNTS_OUTPUT_FOLDER>
					--spec <GENES_OF_INTEREST(comma separated list)>
					"""

#input folder contains kallisto count tables

import os, glob, sys
import matplotlib.pyplot as plt

# --- end of imoprts --- #

def get_samples_per_assembly_ID( species_info_file ):
	"""! @brief get all samples assigned to a certain ID """
	
	samples_per_ID = {}
	with open( species_info_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				if "," in parts[-1]:
					samples_per_ID[ parts[1] ] += parts[-1].split(',')
				else:
					samples_per_ID[ parts[1] ] += [ parts[-1] ]
			except KeyError:
				if "," in parts[-1]:
					samples_per_ID.update( { parts[1]: parts[-1].split(',') } )
				else:
					samples_per_ID.update( { parts[1]: [ parts[-1] ] } )
			line = f.readline()
	return samples_per_ID


def load_expression( data_input_folder, samples ):
	"""! @brief load expression values of all samples """
	
	# --- load data --- #
	final_samples = []
	tpms = {}
	counts = {}
	for sample in samples:
		filename = data_input_folder + sample + ".tsv"
		if os.path.isfile( filename ):
			final_samples.append( sample )
			with open( filename, "r" ) as f:
				f.readline()	#remove header
				line = f.readline()
				while line:
					parts = line.strip().split('\t')
					try:
						counts[ parts[0] ].append( parts[3] )
						tpms[ parts[0] ].append( parts[4] )
					except KeyError:
						counts.update( { parts[0] : [ parts[3] ] } )
						tpms.update( { parts[0] : [ parts[4] ] } )
					line =f.readline()
	
	# --- check tpm and count lengths of all genes --- #
	all_values = tpms.values() + counts.values()
	lengths = []
	for val in all_values:
		lengths.append( len( val ) )
	
	if len( list( set( lengths ) ) ) == 1:
		return tpms, counts, final_samples
	elif len(  list( set( lengths ) ) ) > 1:
		print "ERROR: inconsistent count table among " + str( samples )
		return tpms, counts, []
	else:
		return tpms, counts, []


def main( arguments ):
	"""! @brief run everything """
	
	data_input_folder = arguments[ arguments.index( '--in' )+1 ]
	tpm_folder = arguments[ arguments.index( '--tpms' )+1 ]
	counts_folder = arguments[ arguments.index( '--counts' )+1 ]
	species_info_file = arguments[ arguments.index( '--spec' )+1 ]
	
	if not os.path.exists( counts_folder ):
		os.makedirs( counts_folder )
	if not os.path.exists( tpm_folder ):
		os.makedirs( tpm_folder )
	
	samples_per_ID = get_samples_per_assembly_ID( species_info_file )	#dictionary: assemblyID: [ sample1, sample2, ...]
	for ID in samples_per_ID.keys():
		tpms, counts, samples = load_expression( data_input_folder, samples_per_ID[ ID ] )
		tpm_file = tpm_folder + ID + ".txt"
		count_file = counts_folder + ID + ".txt"
		
		if len( samples ) > 0:
			with open( tpm_file, "w" ) as out:
				out.write( "\t" + "\t".join( samples ) + '\n' )
				for gene in sorted( tpms.keys() ):
					out.write( "\t".join( [ gene ] + tpms[ gene ] ) + '\n' )
			
			with open( count_file, "w" ) as out:
				out.write( "\t" + "\t".join( samples ) + '\n' )
				for gene in sorted( counts.keys() ):
					out.write( "\t".join( [ gene ] + counts[ gene ] ) + '\n' )



if '--in' in sys.argv and '--tpms' in sys.argv and '--counts' in sys.argv and '--spec' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
