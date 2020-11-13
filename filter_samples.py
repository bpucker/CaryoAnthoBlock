### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
						python filter_samples.py
						--in <TPM_FILE>| --indir <TPM_FILE_DIRECTORY>
						--out <OUTPUT_FILE>|--outdir <TPM_OUTPUT_FILE_DIRECTORY>
						--min <MIN_PERCENT_EXPRESSION_ON_TOP100>[10]
						--max <MAX_PERCENT_EXPRESSION_ON_TOP100>[80]
						"""

import os, sys, glob
import matplotlib.pyplot as plt

# --- end of imports --- #


def load_all_TPMs( exp_file ):
	"""! @brief load all values from given TPM file """
	
	data = {}
	genes = []
	with open( exp_file, "r" ) as f:
		headers = f.readline().strip()
		if "\t" in headers:
			headers = headers.split('\t')
		else:
			headers = [ headers ]
		if headers[0] == "gene":
			headers = headers[1:]
		for header in headers:
			data.update( { header: [] } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			genes.append( parts[0] )
			for idx, val in enumerate( parts[1:] ):
				data[ headers[ idx ] ].append( float( val ) )
			line = f.readline()
	return data, genes


def main( arguments ):
	"""! @brief run everything """
	
	try:	#analysing just one file
		input_files = [ arguments[ arguments.index('--in')+1 ] ]
		output_files = [ arguments[ arguments.index('--out')+1 ] ]
		doc_files = [ output_files[0] + ".doc" ]
		fig_files = [ output_files[0] + ".pdf" ]
	except ValueError:	#analysing all files in a folder
		input_folder = arguments[ arguments.index('--indir')+1 ]
		output_folder = arguments[ arguments.index('--outdir')+1 ]
		if not os.path.exists( output_folder ):
			os.makedirs( output_folder )
		if output_folder[-1] != "/":
			output_folder += "/"
		input_files = glob.glob( input_folder + "*.txt" )
		output_files = []
		doc_files = []
		fig_files = []
		for filename in input_files:
			ID = filename.split('/')[-1].split('.')[0]
			output_files.append( output_folder + ID + ".txt" )
			doc_files.append( output_folder + ID + ".txt.doc" )
			fig_files.append( output_folder + ID + ".txt.pdf" )
	
	if '--min' in arguments:
		min_cutoff = int( arguments[ arguments.index('--min')+1 ] )
	else:
		min_cutoff = 10	#value in percent
	if '--max' in arguments:
		max_cutoff = int( arguments[ arguments.index('--max')+1 ] )
	else:
		max_cutoff = 80
	
	# --- run analysis of all data in folder/file --- #
	for k, input_file in enumerate( input_files ):
		output_file = output_files[ k ]
		doc_file = doc_files[ k ]
		fig_file = fig_files[ k ]
		
		valid_samples = []
		with open( doc_file, "w" ) as out:
			out.write( "SampleName\tPercentageOfTop100\tPercentageOfTop500\tPercentageOfTop1000\n" )
			TPM_data, genes = load_all_TPMs( input_file )
			for key in sorted( TPM_data.keys() ):
				new_line = [ key ]
				selection = sorted( TPM_data[ key ] )
				try:
					val = 100.0 * sum( selection[-100:] ) / sum( selection )
				except ZeroDivisionError:
					val = 0
				new_line.append( val )
				if min_cutoff < val < max_cutoff:
					valid_samples.append( key )
				if len( selection ) > 500 and val > 0:
					new_line.append( 100.0 * sum( selection[-500:] ) / sum( selection ) )
				else:
					new_line.append( "n/a" )
				if len( selection ) > 1000 and val > 0:
					new_line.append( 100.0 * sum( selection[-1000:] ) / sum( selection ) )
				else:
					new_line.append( "n/a" )
				out.write( "\t".join( map( str, new_line ) ) + "\n" )
		
		print "number of valid sample: " + str( len( valid_samples ) )
		print "number of invalid sample: " + str( len( TPM_data.keys() ) - len( valid_samples ) )
		
		# --- generate output file --- #
		if len( valid_samples ) > 0:
			with open( output_file, "w" ) as out:
				out.write( "gene\t" + "\t".join( valid_samples )+ "\n" )
				for idx, gene in enumerate( genes ):
					new_line = [ gene ]
					for sample in valid_samples:
						new_line.append( TPM_data[ sample ][ idx ] )
					out.write( "\t".join( map( str, new_line ) ) + "\n" )
		else:
			print "WARNING: no valid samples in data set!"
		# --- generate figure --- #
		values = []
		with open( doc_file, "r" ) as f:
			f.readline()	#remove header
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				values.append( float( parts[1] ) )
				line = f.readline()
		
		values = [ x for x in values if str(x) != 'nan' ]
		
		fig, ax = plt.subplots()
		
		ax.hist( values, bins=100, color="green" )
		ax.set_xlabel( "Percentage of expression on top100 genes" )
		ax.set_ylabel( "Number of analyzed samples" )
		
		fig.savefig( fig_file )
		plt.close( "all" )


if '--in' in sys.argv and '--out' in sys.argv:	#analyse one file
	main( sys.argv )
elif '--indir' in sys.argv and '--outdir' in sys.argv:	#analyse all files in one folder
	main( sys.argv )
else:
	sys.exit( __usage__ )
