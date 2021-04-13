### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python jcvi_wrapper_genes.py
					--gff <COMMA_SEPARATED_LIST_OF_GFF_FILES>
					--cds <COMMA_SEPARATED_LIST_OF_CDS_FILES>
					--feature <COMMA_SEPARATED_LIST_OF_FEATURE_NAMES>(mRNA|transcript)
					--ID <COMMA_SEPARATED_LIST_OF_ID_NAMES>(ID)
					--specs <COMMA_SEPARATED_LIST_OF_SPECIES_NAMES>
					--chromosomes <COMMA_SEPARATED_LIST_OF_CHROMOSOME_NAMES>
					"""

import os, sys, re

# --- end of imports --- #

def prepare_customized_blocks_file( genes, input_file, output_file ):
	"""! @brief select blocks of interest """
	
	data = []
	with open( input_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				#ID = re.findall( "[a-z]{4}", parts[0] )[0]
				data.append( { 'ID': parts[0], 'line': line } )
			except IndexError:
				print line
			line = f.readline()
	
	# --- get min and max index --- #
	min_index = len( data )
	max_index = 0
	for gene in genes:
		for idx, entry in enumerate( data ):
			if gene == entry['ID']:
				if idx < min_index:
					min_index = idx + 0
				if idx > max_index:
					max_index = idx + 0
	
	if min_index < max_index:
		with open( output_file, "w" ) as out:
			lines = []
			for each in data[ min_index: max_index+1 ]:
				lines.append( each['line'] )
			out.write( "".join( lines ) )
	else:
		print "ERROR (BP): min index >= max index (candidate genes missing?)"
	

def merge_block_files( block_file_per_spec, merged_block_file ):
	"""! @brief merge single block files """
	
	# --- read all data --- #
	mapping_table = {}
	sample_names = []
	genes = []
	for idx, filename in enumerate( block_file_per_spec ):
		sample_names.append( filename )
		tmp = {}
		with open( filename, "r" ) as f:
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				if idx == 0:
					genes.append( parts[0] )
				tmp.update( { parts[0]: parts[1:] } )
				line = f.readline()
		mapping_table.update( { filename: tmp } )
	# --- generate output file --- #
	with open( merged_block_file, "w" ) as out:
		for gene in genes:
			new_line = [ gene ]
			for sample in sample_names:
				new_line.append( "\t".join( mapping_table[ sample ][ gene ] ) )
			out.write( "\t".join( new_line ) + "\n" )


def modify_lifted_anchor_file( lifted_anchor_file ):
	"""! @brief modify lifted anchor file to remove entries of empty contigs """
	
	with open( lifted_anchor_file, "r" ) as f:
		content = f.read()
	content = content.replace( "###\n###\n", "###\n" ).replace( "###\n###\n", "###\n" ).replace( "###\n###\n", "###\n" ).replace( "###\n###\n", "###\n" ).replace( "###\n###\n", "###\n" ).replace( "###\n###\n", "###\n" ).strip()
	
	with open( lifted_anchor_file, "w" ) as out:
		out.write( content )


def main( arguments ):
	"""! @brief run everything """
	
	gff_files = arguments[ arguments.index('--gff')+1 ].split(',')
	cds_files = arguments[ arguments.index('--cds')+1 ].split(',')
	feature_types = arguments[ arguments.index('--feature')+1 ].split(',')	#tag in third column
	ID_types = arguments[ arguments.index('--ID')+1 ].split(',')	#tag in last field
	spec_names = arguments[ arguments.index('--specs')+1 ].split(',') #species names appear as labels in figure
	chromosomes_of_interst = []
	raw_chromosomes = arguments[ arguments.index('--chromosomes')+1 ].split(',')
	for each in raw_chromosomes:
		if "@@" in each:	#multiple chromosomes per species are possible (?)
			chromosomes_of_interst.append( each.split('@@') )
		else:
			chromosomes_of_interst.append( [ each ] )

	geneIDs = arguments[ arguments.index('--genes')+1 ].split(',')	#BeetSet2 IDs

	cscore_cutoff = 0.1	#needs to be a value between 0 and 1
	number_of_matching_regions = 1	#get one corresponding region for each region in reference (needs to be increased)

	output_folder = arguments[ arguments.index('--out')+1 ]

	# --- generate output folder and set it as working directory (important, because relative paths are used!) --- #
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )

	os.chdir( output_folder )


	# --- production of bed files --- #
	bed_files = []
	for idx, filename in enumerate( gff_files ):
		bed_file = spec_names[ idx ] + ".bed"
		if not os.path.isfile( bed_file ):
			os.popen( "python -m jcvi.formats.gff bed --type " + feature_types[ idx ] + " --key=" + ID_types[ idx ] + " " + filename + " -o " + bed_file )
		bed_files.append( bed_file )


	# --- generate clean CDS files --- #
	clean_cds_files = []
	for idx, filename in enumerate( cds_files ):
		cds_file = spec_names[ idx ] + ".cds"
		if not os.path.isfile( cds_file ):
			os.popen( "python -m jcvi.formats.fasta format " + filename + " " + cds_file )
		clean_cds_files.append( cds_file )


	# --- generate one block file per species --- #
	block_file_per_spec = []
	for idx, spec in enumerate( spec_names ):
		if idx > 0:
			# --- identification of putative orthologs --- #
			anchor_file = ".".join( [ spec_names[0], spec ] ) + ".anchors "
			lifted_anchor_file = ".".join( [ spec_names[0], spec ] ) + ".lifted.anchors"
			if not os.path.isfile( lifted_anchor_file ):
				os.popen( "python -m jcvi.compara.catalog ortholog " + " ".join( [ spec_names[0], spec ] ) + " --cscore=" + str( cscore_cutoff ) + " --no_strip_names" )

			# --- generate more succint anchors file --- #
			modify_lifted_anchor_file( lifted_anchor_file )
			block_file = ".".join( [ spec_names[0], spec ] ) + ".i" + str( number_of_matching_regions ) + ".blocks"
			if not os.path.isfile( block_file ):
				os.popen( "python -m jcvi.compara.synteny mcscan " + bed_files[0] + " " + lifted_anchor_file + " --iter=" + str( number_of_matching_regions ) + " -o " + block_file )
			block_file_per_spec.append( block_file )

	# --- merge block files of all species --- #
	merged_block_file = spec_names[0] + ".blocks"
	#os.popen( "python -m jcvi.formats.base join " + " ".join( block_file_per_spec ) + " --noheader | cut -f1,2,4,6 > " + merged_block_file )
	merge_block_files( block_file_per_spec, merged_block_file )	#MY SOLUTION, but might enfores plot with wrong locus

	# --- just select a random region --- #
	#os.popen( "head -50 " + block_file + " > blocks" )
	prepare_customized_blocks_file( geneIDs, merged_block_file, "blocks" )


	# --- generate seqids and layout file --- #
	layout_file = "blocks.layout"

	yvalues = [ "0.6", "0.8",  "0.8", " .4", " .4" ]
	xvalues = [ "0.5", "0.3", "0.7", "0.3", "0.7" ]
	has = [ "center", "center", "center", "center", "center" ]
	vas = [ "bottom", "bottom", "bottom", "bottom", "bottom" ]
	with open( layout_file, "w" ) as out:
		out.write( "# x, y, rotation, ha, va, color, ratio, label\n" )
		for idx, spec in enumerate( spec_names ):
			out.write( xvalues[ idx ] + ", " + yvalues[ idx ] + ", 0, " + has[ idx ] + ", " + vas[idx] + ", m,  1," + spec + ", " + chromosomes_of_interst[ idx ][0] + "\n" )
		out.write( "# edges\n" )
		for idx, spec in enumerate( spec_names[1:] ):
			out.write( "e, 0, " + str( idx+1 ) + "\n" )

	merged_bed_file = ".".join( spec_names ) + ".bed"
	os.popen( "cat " + " ".join( bed_files ) + " > " + merged_bed_file )

	os.popen( "python -m jcvi.graphics.synteny blocks " + merged_bed_file + " " + layout_file + " --glyphstyle=arrow --genelabelsize 4 --scalebar" )


	print "COMMAND TO RE-RUN PLOT AFTER ADJUSTING LAYOUT FILE:\n\n"
	print "python -m jcvi.graphics.synteny blocks " + merged_bed_file + " " + layout_file + " --glyphstyle=arrow --genelabelsize 4 --scalebar"
	print "\n\n"


if '--gff' in sys.argv and '--cds' in sys.argv and '--feature' in sys.argv and '--ID' in sys.argv and '--specs' in sys.argv and '--chromosomes' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
