### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###


__usage__ = """
					python black_list_cleaning.py
					--cdsin <CDS_INPUT_FOLDER>
					--pepin <PEP_INPUT_FOLDER>
					--black <BLACK_LIST>
					--cdsout <CDS_OUTPUT_FOLDER>
					--pepout <PEP_OUTPUT_FOLDER>
					"""

import glob, os, sys

# --- end of imports --- #

def load_black_list( black_list_file ):
	"""! @brief load all sequence IDs from given file (one per line) """
	
	black_list = {}
	with open( black_list_file, "r" ) as f:
		line = f.readline()
		while line:
			black_list.update( { line.strip(): None } )
			line = f.readline()
	return black_list


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def clean_files( input_folder, output_folder, black ):
	"""! @brief clean all files in provided folder """
	
	filenames = glob.glob( input_folder + "*.fa" ) + glob.glob( input_folder + "*.fasta" )
	for filename in filenames:
		output_file = output_folder + filename.split('/')[-1]
		seqs = load_sequences( filename )
		with open( output_file, "w" ) as out:
			for key in seqs.keys():
				try:
					black[ key ]
				except KeyError:
					out.write( '>' + key + "\n" + seqs[ key ] + "\n" )


def main( arguments ):
	"""! @brief run everything """

	cds_input_folder = arguments[ arguments.index('--cdsin')+1 ]
	pep_input_folder = arguments[ arguments.index('--pepin')+1 ]
	black_list_file = arguments[ arguments.index('--black')+1 ]

	cds_output_folder = arguments[ arguments.index('--cdsout')+1 ]
	pep_output_folder = arguments[ arguments.index('--pepout')+1 ]
	
	if not os.path.exists( cds_output_folder ):
		os.makedirs( cds_output_folder )
	
	if not os.path.exists( pep_output_folder ):
		os.makedirs( pep_output_folder )
	
	
	black_list = load_black_list( black_list_file )
	
	clean_files( cds_input_folder, cds_output_folder, black_list )
	clean_files( pep_input_folder, pep_output_folder, black_list )


if "--cdsin" in sys.argv and '--pepin' in sys.argv and '--black' in sys.argv and '--cdsout' in sys.argv and '--pepout' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
