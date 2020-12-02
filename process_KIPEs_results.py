### Boas Pucker ###
### bp423@cam.ac.uk ###
### v0.13 ###

__usage__ = """
					python process_KIPEs_results.py
					--in <LIST_OF_INTPUT_FOLDERS(comma separated list)>
					--out <OUTPUT_FOLDER>
					--genes <GENES_OF_INTEREST(comma separated list)>
					--roi <UNDERLINE_SEP_GENES_COMMA_SEP_RESIDUES_OF_INTEREST>[False]
					"""

import os, glob, sys
import matplotlib.pyplot as plt

# --- end of imoprts --- #


def get_info_per_species( inputfile ):
	"""! @brief load info per sequence """
	
	infos = {}
	with open( inputfile, "r" ) as f:
		checked_residues = f.readline().strip().split('\t')[1:]	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			tmp = {}
			for idx, part in enumerate( parts[1:] ):
				tmp.update( { checked_residues[ idx ]: part } )
			infos.update( { parts[0]: tmp } )
			line = f.readline()
	return infos, checked_residues


def generate_aa_cons_fig( summary_file, fig_file, gene, residues_of_interest ):
	"""! @brief generate barplot to illustrate conservation of different amino acids across species """
	
	# --- load data from file --- #
	with open( summary_file, "r" ) as f:
		cons_res = f.readline().strip().split('\t')[1:]	#header
		counter = {}
		for aa in cons_res:
			counter.update( { aa: 0 } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			for idx, part in enumerate( parts[1:] ):
				if part == "True":
					counter[ cons_res[ idx ] ] += 1
			line = f.readline()
	
	if not residues_of_interest:
		residues_of_interest = cons_res
	
	# --- generate plot --- #
	x_values = []
	y_values = []
	for aa in cons_res:
		if aa in residues_of_interest:
			x_values.append( len( x_values ) )
			y_values.append( counter[ aa ] )
	
	fig, ax = plt.subplots( figsize=(15, 5) )
	
	ax.bar( x_values, y_values, tick_label=residues_of_interest )
	
	ax.set_ylabel( "number of sequences" )
	ax.set_xlabel( "amino acid position" )
	
	ax.set_title( gene )
	plt.subplots_adjust( left=0.05, right=0.999, top=0.95 )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	
	fig.savefig( fig_file )
	plt.close( "all" )


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


def main( arguments ):
	"""! @brief run everything """
	
	input_folders = arguments[ arguments.index( '--in' )+1 ]
	if "," in input_folders:
		input_folders = input_folders.split(',')
	else:
		input_folders = [ input_folders ]
	genes_of_interest = arguments[ arguments.index( '--genes' )+1 ]
	if "," in genes_of_interest:
		genes_of_interest = genes_of_interest.split(',')
	else:
		genes_of_interest = [ genes_of_interest ]
	
	output_folder = arguments[ arguments.index( '--out' )+1 ]
	
	if '--roi' in arguments:
		tmp_all_residues_of_interest = arguments[ arguments.index( '--roi' )+1 ]
		print tmp_all_residues_of_interest
		if "_" in tmp_all_residues_of_interest:
			gene_parts = tmp_all_residues_of_interest.split('_')
		else:
			gene_parts = [ tmp_all_residues_of_interest ]
		all_residues_of_interest = []
		for each in gene_parts:
			if ',' in each:
				all_residues_of_interest.append( each.split(',') )
			else:
				all_residues_of_interest.append( [ each ] )
	else:
		all_residues_of_interest = [ False ] * len( genes_of_interest )
	
	print "number of residues to check: " + str( len( all_residues_of_interest[0] ) )
	
	
	if '--tolerate' in arguments:
		tolerate = int( arguments[ arguments.index( '--tolerate' )+1 ] )
	else:
		tolerate = 0

	if '--gaps' in arguments:
		gaps_allowed = True
	else:
		gaps_allowed = False

	if output_folder[-1] != '/':
		output_folder += "/"

	seq_output_folder = output_folder + "seq/"
	info_output_folder = output_folder + "info/"
	cons_output_folder = output_folder + "cons/"
	
	if not os.path.exists( seq_output_folder ):
		os.makedirs( seq_output_folder )
	if not os.path.exists( info_output_folder ):
		os.makedirs( info_output_folder )
	if not os.path.exists( cons_output_folder ):
		os.makedirs( cons_output_folder )
	
	# --- generating data summary files --- #
	seq_output_files = {}
	info_output_files = {}
	seq_output_file_handles = {}
	info_output_file_handles = {}
	for gene in genes_of_interest:
		seq_output_files.update( { gene: seq_output_folder + gene + ".fasta" } )
		info_output_files.update( { gene: info_output_folder + gene + ".txt" } )
		seq_output_file_handles.update( { gene: open( seq_output_folder + gene + ".fasta", "w" ) } )
		info_output_file_handles.update( { gene: open( info_output_folder + gene + ".txt", "w" ) } )
	
	pass_seq_counter = 0
	header_status_per_gene = {}
	for gene in genes_of_interest:
		header_status_per_gene.update( { gene: True } )
	for idx1, folder in enumerate( input_folders ):
		sub_folders = sorted( glob.glob( folder + "*/" ) )
		for idx2, subf in enumerate( sub_folders ):
			for g, gene in enumerate( genes_of_interest ):
				info_file_name = subf + "conserved_residues/" + gene + "_conserved_residues.txt"
				fasta_file = subf + "final_pep_files/" + gene + ".fasta"
				if os.path.isfile( info_file_name ):
					info, residues = get_info_per_species( info_file_name )
					seqs = load_sequences( fasta_file )
					if header_status_per_gene[ gene ]:
						info_output_file_handles[ gene ].write( "sequences\t" + "\t".join( residues ) + '\n' )
						header_status_per_gene[ gene ] = False
					for seq_info in info.keys():
						new_line = [ seq_info ]
						for aa in residues:
							new_line.append( info[ seq_info ][ aa ] )
						info_output_file_handles[ gene ].write( "\t".join( new_line ) + '\n' )
						ok_counter = 0
						for roi in all_residues_of_interest[ g ]:
							if info[ seq_info ][ roi ] == "True":
								ok_counter += 1
							elif gaps_allowed:
								if info[ seq_info ][ roi ] == "-":
									ok_counter += 1
						#if "True" in new_line[1:] and len( list( set( new_line[1:]  ) ) ) == 1:
						if ok_counter >= len( all_residues_of_interest[ g ] )-tolerate:
							seq_output_file_handles[ gene ].write( '>' + new_line[0] + '\n' + seqs[ new_line[0] ] + '\n' )
							pass_seq_counter += 1
	print "number of passed genes: " + str( pass_seq_counter )
	
	# --- running analysis on data --- #
	for x, gene in enumerate( genes_of_interest ):
		summary_file = info_output_files[ gene ]
		fig_file = cons_output_folder + gene + ".pdf"
		generate_aa_cons_fig( summary_file, fig_file, gene, all_residues_of_interest[ x ] )
		residues_of_interest = False


if '--in' in sys.argv and '--out' in sys.argv and '--genes' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
