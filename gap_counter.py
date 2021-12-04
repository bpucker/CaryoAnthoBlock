### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.1 ###

__usage__ = """
					python3 gap_counter.py
					--in <INPUT_FILE>
					bug reports and feature requests: b.pucker@tu-braunschweig.de
					"""

import os, sys, re

# --- end of imports --- #


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(" ")[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(" ")[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	
	seqs = load_sequences( input_file )
	
	for key in list( seqs.keys() ):
		seq = seqs[ key ].upper()
		gaps = re.findall( "[N]+", seq )
		sys.stdout.write( key + ": " + str( len( gaps ) ) + "\n" )
		sys.stdout.flush()


if '--in' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
