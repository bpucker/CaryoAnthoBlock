### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python kallisto_pipeline.py
					--cds <CDS_FILE_FOLDER>
					--info <SPECIES_AND_DATA_TABLE>
					--reads <FOLDER_CONTAINING_FASTQs>
					--tmp <TMP_FOLDER>
					--out <OUTPUT_FOLDER>
					
					optional:
					--kallisto <FULL_PATH_TO_KALLISTO>[kallisto]
					--cpus <NUMBER_OF_THREADS>[10]
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""


import os, sys, glob


# --- end of imports --- #


def load_data( input_file ):
	"""! @brief extract all IDs from given input file """
	
	all_data = []
	with open( input_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		index = 0
		while line:
			parts = line.strip().split('\t')
			all_data.append( { 'project': parts[0], 'ID': parts[1], 'taxon': parts[2], 'status': parts[3], 'SRA': parts[4].split(','), 'idx': index } )
			index += 1
			line = f.readline()
	print "number of species: " + str( len( all_data ) )
	return all_data


def construct_mapping_table( data ):
	"""! @brief construct mapping table for all data """
	
	mapping_table = {}
	for entry in data:
		for ID in entry['SRA']:
			mapping_table.update( { ID: entry } )
	return mapping_table


def get_data_for_jobs_to_run( read_file_folders, data_per_ID, final_output_folder, cds_file_input_folder, tmp_cluster_folder ):
	"""! @brief collect all infos to run jobs """
	
	jobs_to_do = []
	for folder in read_file_folders:
		ID = folder.split('/')[-1]
		status = True
		
		# --- get read file --- #
		PE_status = True
		read_file1 = folder + "/" + ID + "_pass_1.fastq.gz"
		if not os.path.isfile( read_file1 ):
			print "ERROR: file missing - " + read_file1
			PE_status = False
		read_file2 = folder + "/" + ID + "_pass_2.fastq.gz"
		if not os.path.isfile( read_file2 ):
			print "ERROR: file missing - " + read_file2
			PE_status = False
		if not PE_status:
			read_file1 = folder + "/" + ID + "_pass.fastq.gz"
			if not os.path.isfile( read_file1 ):
				status = False
			read_file2 = False
		
		# --- get reference for quantification --- #
		cds_file = cds_file_input_folder + data_per_ID[ ID ]['ID'] + ".cds.fa"
		if not os.path.isfile( cds_file ):
			print "ERROR: file missing - " + cds_file
			status = False
		output_dir = tmp_cluster_folder + ID + "/"
		if not os.path.exists( output_dir ):
			os.makedirs( output_dir )
		index_file = output_dir + "index"
		tmp_result_file = output_dir + "abundance.tsv"
		final_result_file = final_output_folder + ID + ".tsv"
		if os.path.isfile( final_result_file ):
			status = False
		if status:
			jobs_to_do.append( { 'r1': read_file1, 'r2': read_file2, 'cds': cds_file, 'out': output_dir, 'index': index_file, 'tmp': tmp_result_file, 'fin': final_result_file, "ID": ID } )
	return jobs_to_do


def job_executer( jobs_to_run, kallisto, threads ):
	"""! @brief run all jobs in list """
	
	for idx, job in enumerate( jobs_to_run ):
		print "running job " + str( idx+1 ) + "/" + str( len( jobs_to_run ) ) + " - " + job["ID"]
		cmd1 = " ".join( [ kallisto, "index", "--index="+job['index'], "--make-unique", job["cds"] ] )
		os.popen( cmd1 )
		
		if job['r2']:
			cmd2 = " ".join( [ kallisto, "quant", "--index="+job['index'], "--output-dir="+job['out'], "--threads "+str( threads ), job['r1'], job['r2'] ] )
		else:
			cmd2 = " ".join( [ kallisto, "quant", "--index="+job['index'], "--single -l 600 -s 300", "--output-dir="+job['out'], "--threads "+str( threads ), job['r1'] ] )
		os.popen( cmd2 )
		
		os.popen( "cp " + job["tmp"] + " " + job["fin"] )


def main( arguments ):
	"""! @brief run everything """

	cds_file_input_folder = arguments[ arguments.index('--cds')+1 ]
	data_table = arguments[ arguments.index('--info')+1 ]
	read_file_folders = arguments[ arguments.index('--reads')+1 ]
	if "," in read_file_folders:
		read_file_folders = read_file_folders.split(',')
	else:
		read_file_folders = [ read_file_folders ]
	tmp_cluster_folder = arguments[ arguments.index('--tmp')+1 ]
	final_output_folder = arguments[ arguments.index('--out')+1 ]
	
	if '--kallisto' in arguments:
		kallisto = arguments[ arguments.index('--kallisto')+1 ]
	else:
		kallisto = "kallisto"
	if '--cpus' in arguments:
		threads = int( arguments[ arguments.index('--cpus')+1 ] )
	else:
		threads = 10
	
	# --- load data --- #
	data = load_data( data_table )
	data_per_ID = construct_mapping_table( data )	#mapping table of IDs to all data info of that sample
	single_read_file_folders = []
	for read_file_folder in read_file_folders:
		single_read_file_folders += [ x[0] for x in os.walk( read_file_folder ) ][1:]
	print "number of FASTQ file folders detected: " + str( len( single_read_file_folders ) )
	
	# --- prepare jobs to run --- #
	jobs_to_run = get_data_for_jobs_to_run( single_read_file_folders, data_per_ID, final_output_folder, cds_file_input_folder, tmp_cluster_folder )
	print "number of jobs to run: " + str( len( jobs_to_run ) )
	
	# --- run jobs --- #
	job_executer( jobs_to_run, kallisto, threads )


if '--cds' in sys.argv and '--info' in sys.argv and '--reads' in sys.argv and '--tmp' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
