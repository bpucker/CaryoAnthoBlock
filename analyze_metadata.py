### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.1 ###

__usage__ = """
					python3 analyze_metadata.py
					--sra <SRA_META_DATA_FILE>
					--summary <CARYO_SUMMARY_FILE>
					--out <OUTPUT_FOLDER>
					bug reports and feature requests: b.pucker@tu-braunschweig.de
					"""

import os, sys, re
from operator import itemgetter

# --- end of imports --- #

def load_study_per_run( sra_file ):
	"""! @brief load study ID per run ID """
	
	run2study = {}
	with open( sra_file, "r" ) as f:
		headers = f.readline().strip().split(',')
		line = f.readline()
		while line:
			parts = line.strip().split(',')
			#run2study.update( { parts[0]: parts[ headers.index( "SRA Study" ) ] } )
			try:
				run2study.update( { parts[0]: re.findall( "SRP\d+", line )[0] } )
			except IndexError:
				try:
					run2study.update( { parts[0]: re.findall( "ERP\d+", line )[0] } )
				except IndexError:
					try:
						run2study.update( { parts[0]: re.findall( "DRP\d+", line )[0] } )
					except IndexError:
						print( line[:100] )
			line = f.readline()
	return run2study


def load_sample_per_run( sra_file ):
	"""! @brief load study ID per run ID """
	
	run2anno = {}
	with open( sra_file, "r" ) as f:
		headers = f.readline().strip().split(',')
		line = f.readline()
		while line:
			parts = line.strip().split(',')
			#run2anno.update( { parts[0]: parts[ headers.index( "tissue" ) ] } )	#tissue
			line = line.lower()
			pot_tissue = []
			if "leaf" in line:
				pot_tissue.append( "leaf" )
			elif "leave" in line:	#to catch plural (and typos)
				pot_tissue.append( "leaf" )
			if "flower" in line:
				pot_tissue.append( "flower" )
			if "stem" in line:
				pot_tissue.append( "stem" )
			if "root" in line:
				pot_tissue.append( "root" )
			if "seedling" in line:
				pot_tissue.append( "seedling" )
			elif "seed" in line:
				pot_tissue.append( "seed" )
			
			if len( pot_tissue ) == 1:
				run2anno.update( { parts[0]: pot_tissue[0] } )
			else:
				run2anno.update( { parts[0]: "unknown" } )	#optional (could also be ignored)
			line = f.readline()
	return run2anno


def load_taxon_file( summary_file ):
	"""! @brief load information from taxon file """
	
	runs_per_species, pigmentation_per_species, ID2spec = {}, {}, {}
	with open( summary_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 4:
				if len( parts[4] ) > 3:
					if "," in parts[4]:
						runs = parts[4].split(',')
					else:
						runs = [ parts[4] ]
					runs_per_species.update( { parts[0]: runs } )
					pigmentation_per_species.update( { parts[0]: parts[3] } )
					ID2spec.update( { parts[0]: parts[2] } )
			line = f.readline()
	return runs_per_species, pigmentation_per_species, ID2spec


def analyze_sample( values, tissue_to_check ):
	"""! @brief analyze distribution of samples """
	
	results = {}
	for tissue in tissue_to_check:
		results.update( { tissue: values.count( tissue ) } )
	return results


def main( arguments ):
	"""! @brief run everything """
	
	sra_file = arguments[ arguments.index('--sra')+1 ]
	summary_file = arguments[ arguments.index('--summary')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- load data --- #
	run2study = load_study_per_run( sra_file )
	run2anno = load_sample_per_run( sra_file )
	runs_per_species, pigmentation_per_species, ID2spec = load_taxon_file( summary_file )
	
	
	# --- studies per species --- #
	number_of_studies_per_species = {}
	for key in list( runs_per_species.keys() ):
		studies = []
		for run in runs_per_species[ key ]:
			try:
				studies.append( run2study[ run ] )
			except KeyError:
				pass
		number_of_studies_per_species.update( { key: list( set( studies ) ) } )
	
	studies_per_species_output_file = output_folder + "studies_per_species.txt"
	with open( studies_per_species_output_file, "w" ) as out:
		for spec in list( number_of_studies_per_species.keys() ):
			if len( number_of_studies_per_species[spec] ) > 0:
				out.write( "\t".join( map( str, [ spec, ID2spec[ spec ], len( number_of_studies_per_species[spec] ), ",".join( number_of_studies_per_species[spec] ) ] ) ) + "\n" )
	
	
	# --- compare tissue distribution between origins at the species level --- #
	most_abundant_tissue_per_species = {}
	for spec in list( runs_per_species.keys() ):
		annos = []
		missing = 0
		for run in runs_per_species[ spec ]:
			try:
				annos.append( run2anno[ run ] )
			except KeyError:
				missing += 1
		#print( spec + ": " + str( len( annos ) ) + "/" + str( missing ) )	#show proportion of annotated samples
		if len( annos ) > 0:
			sample_frequency = []
			for each in list( set( annos ) ):
				sample_frequency.append( { 'anno': each, 'counts': annos.count( each ) } )
			sample_frequency = sorted( sample_frequency, key=itemgetter('counts') )[::-1]
			if len( sample_frequency ) == 1:
				most_abundant_tissue_per_species.update( { spec: sample_frequency[0]['anno'] } )
			else:
				if sample_frequency[0]['counts'] > sample_frequency[1]['counts']:
					most_abundant_tissue_per_species.update( { spec: sample_frequency[0]['anno'] } )
				else:
					pass	#print( annos )
	A, B2, B3, B4 = [], [], [], []
	A_group, B_group = [], []
	for spec in list( pigmentation_per_species.keys() ):
		try:
			tissue = most_abundant_tissue_per_species[ spec ]
			if pigmentation_per_species[ spec ] == "A":
				A.append( tissue )
				A_group.append( tissue )
			elif pigmentation_per_species[ spec ] == "B2":
				B2.append( tissue )
				B_group.append( tissue )
			elif pigmentation_per_species[ spec ] == "B3":
				B3.append( tissue )
				B_group.append( tissue )
			elif pigmentation_per_species[ spec ] == "B4":
				B4.append( tissue )
				B_group.append( tissue )
		except KeyError:
			pass
	
	tissue_to_check = [ "leaf", "flower", "root", "seedling", "stem", "seed", "unknown" ]
	rows = { "all anthocyanins": analyze_sample( A_group, tissue_to_check ),
					"all betalains": analyze_sample( B_group, tissue_to_check ),
					"A": analyze_sample( A, tissue_to_check ),
					"B2": analyze_sample( B2, tissue_to_check ),
					"B3": analyze_sample( B3, tissue_to_check ),
					"B4": analyze_sample( B4, tissue_to_check )
				}
	
	sample_distr_result_file = output_folder + "sample_distr_result_file_species_level.txt"
	with open( sample_distr_result_file, "w" ) as out:
		out.write( "Group\t" + "\t".join( tissue_to_check ) + "\n" )
		for row in list( rows.keys() ):
			new_line = [ row ]
			for tissue in tissue_to_check:
				new_line.append( rows[ row ][ tissue ] )
			out.write( "\t".join( list( map( str, new_line ) ) ) + "\n" )
	
	# --- compare tissue distribution between origins at sample (run) level --- #
	A, B2, B3, B4 = [], [], [], []
	A_group, B_group = [], []
	for spec in list( runs_per_species.keys() ):
		for run in runs_per_species[ spec ]:
			try:
				tissue = run2anno[ run ]
				if pigmentation_per_species[ spec ] == "A":
					A.append( tissue )
					A_group.append( tissue )
				elif pigmentation_per_species[ spec ] == "B2":
					B2.append( tissue )
					B_group.append( tissue )
				elif pigmentation_per_species[ spec ] == "B3":
					B3.append( tissue )
					B_group.append( tissue )
				elif pigmentation_per_species[ spec ] == "B4":
					B4.append( tissue )
					B_group.append( tissue )
			except KeyError:
				pass
	
	rows = { "all anthocyanins": analyze_sample( A_group, tissue_to_check ),
					"all betalains": analyze_sample( B_group, tissue_to_check ),
					"A": analyze_sample( A, tissue_to_check ),
					"B2": analyze_sample( B2, tissue_to_check ),
					"B3": analyze_sample( B3, tissue_to_check ),
					"B4": analyze_sample( B4, tissue_to_check )
				}
	
	sample_distr_result_file = output_folder + "sample_distr_result_file_run_level.txt"
	with open( sample_distr_result_file, "w" ) as out:
		out.write( "Group\t" + "\t".join( tissue_to_check ) + "\n" )
		for row in list( rows.keys() ):
			new_line = [ row ]
			for tissue in tissue_to_check:
				new_line.append( rows[ row ][ tissue ] )
			out.write( "\t".join( list( map( str, new_line ) ) ) + "\n" )


if '--sra' in sys.argv and '--summary' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
