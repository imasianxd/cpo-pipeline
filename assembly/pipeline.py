#!/usr/bin/env python

'''
This script is a wrapper for module one of the cpo-pipeline including QC and Assembly.
It uses Mash2.0, Kraken2.0 and fastqc to check for sequence contamination, quality information and identify a reference genome.
Then attempts to assemble the reads, attempting to filter contamination away if required.

Example usage:

  pipeline.py -i BC11-Kpn005_S2 -f BC11-Kpn005_S2_L001_R1_001.fastq.gz -r BC11-Kpn005_S2_L001_R2_001.fastq.gz -o output_dir -e "Klebsiella pneumoniae" 

Requires pipeline_qc.sh, pipeline_assembly.sh, pipeline_assembly_contaminant.sh. where these scripts are located can be specified with -k. 
'''

import subprocess
import optparse
import os
import datetime
import sys
import time
import urllib.request
import gzip
import collections
import json
import configparser

from parsers import result_parsers

def execute(command, curDir):
	process = subprocess.Popen(command, shell=False, cwd=curDir, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

	# Poll process for new output until finished
	while True:
		nextline = process.stdout.readline()
		if nextline == '' and process.poll() is not None:
			break
		sys.stdout.write(nextline)
		sys.stdout.flush()

	output = process.communicate()[0]
	exitCode = process.returncode

	if (exitCode == 0):
		return output
	else:
		raise subprocess.CalledProcessError(exitCode, command)

def httpGetFile(url, filepath=""):
	if (filepath == ""):
		return urllib.request.urlretrieve(url)
	else:
		urllib.request.urlretrieve(url, filepath)
		return True

def gunzip(inputpath="", outputpath=""):
	if (outputpath == ""):
		with gzip.open(inputpath, 'rb') as f:
			gzContent = f.read()
		return gzContent
	else:
		with gzip.open(inputpath, 'rb') as f:
			gzContent = f.read()
		with open(outputpath, 'wb') as out:
			out.write(gzContent)
		return True

def main():
	
	config = configparser.ConfigParser()
	config.read(os.path.dirname(os.path.realpath(sys.argv[0])) + '/config.ini')
	
	#parses some parameters
	parser = optparse.OptionParser("Usage: %prog [options] arg1 arg2 ...")
	parser.add_option("-i", "--id", dest="id", type="string", help="identifier of the isolate")	
	parser.add_option("-f", "--forward", dest="R1", type="string", help="absolute file path forward read (R1)")
	parser.add_option("-r", "--reverse", dest="R2", type="string", help="absolute file path to reverse read (R2)")
	parser.add_option("-m", "--mash-genomedb", dest="mashGenomeRefDB", default = config['databases']['mash-genomedb'], type="string", help="absolute path to mash reference database")
	parser.add_option("-n", "--mash-plasmiddb", dest="mashPlasmidRefDB", default = config['databases']['mash-plasmiddb'], type="string", help="absolute path to mash reference database")
	parser.add_option("-z", "--kraken2-genomedb", dest="kraken2GenomeRefDB", default = config['databases']['kraken2-genomedb'], type="string", help="absolute path to kraken reference database")
	parser.add_option("-v", "--kraken2-plasmiddb", dest="kraken2PlasmidRefDB", default = config['databases']['kraken2-plasmiddb'], type="string", help="absolute path to kraken reference database")
	parser.add_option("-x", "--busco-db", dest="buscodb", default = config['databases']['busco-db'], type="string", help="absolute path to busco reference database")

	parser.add_option("-o", "--output", dest="output", default='./', type="string", help="absolute path to output folder")
	parser.add_option("-k", "--script-path", dest="scriptDir", default=config['scripts']['script-path'], type="string", help="absolute file path to this script folder")

	#used for parsing 
	parser.add_option("-e", "--expected", dest="expectedSpecies", default="NA/NA/NA", type="string", help="expected species of the isolate")
	
	#parallelization, useless, these are hard coded to 8cores/64G RAM
	#parser.add_option("-t", "--threads", dest="threads", default=8, type="int", help="number of cpu to use")
	#parser.add_option("-p", "--memory", dest="memory", default=64, type="int", help="memory to use in GB")

	(options,args) = parser.parse_args()

	curDir = os.getcwd()
	outputDir = options.output
	mashdb = options.mashGenomeRefDB
	mashplasmiddb=options.mashPlasmidRefDB
	kraken2db = options.kraken2GenomeRefDB
	kraken2plasmiddb=options.kraken2PlasmidRefDB
	expectedSpecies = options.expectedSpecies
	#threads = options.threads
	#memory = options.memory
	tempDir = outputDir + "/shovillTemp"
	scriptDir = options.scriptDir
	buscodb = options.buscodb
	ID = options.id
	R1 = options.R1
	R2 = options.R2

	
	notes = []
	#init the output list
	output = []
	
	print(str(datetime.datetime.now()) + "\n\nID: " + ID + "\nR1: " + R1 + "\nR2: " + R2)
	output.append(str(datetime.datetime.now()) + "\n\nID: " + ID + "\nR1: " + R1 + "\nR2: " + R2)

	print("step 1: preassembly QC")

	print("running pipeline_qc.sh")
	#input parameters: 1 = id, 2= forward, 3 = reverse, 4 = output, 5=mashgenomerefdb, $6=mashplasmidrefdb, $7=kraken2db, $8=kraken2plasmiddb
	cmd = [scriptDir + "/pipeline_qc.sh", ID, R1, R2, outputDir, mashdb, mashplasmiddb, kraken2db, kraken2plasmiddb]
	result = execute(cmd, curDir)

	print("Parsing the QC results")
	
	#parse read stats
	'''
	pathToMashLog = outputDir + "/qcResult/" + ID + "/" + "mash.log"
	pathToTotalBP = outputDir + "/qcResult/" + ID + "/" + "totalbp"
	stats = result_parsers.parse_read_stats(pathToMashLog, pathToTotalBP)
	'''
	
	#parse genome mash results
	pathToMashGenomeScreenTSV = outputDir + "/qcResult/" + ID + "/" + "mashscreen.genome.tsv"
	mash_hits = result_parsers.parse_mash_result(pathToMashGenomeScreenTSV)

	# 'shared_hashes' field is string in format '935/1000'
	# Threshold is 300 below highest numerator (ie. '935/100' -> 635)
	mash_hits_score_threshold = int(mash_hits[0]['shared_hashes'].split("/")[0]) - 300
	print("*** mash_hits_score_threshold: " + str(mash_hits_score_threshold))
	def score_above_threshold(mash_result, score_threshold):
		score = int(mash_result['shared_hashes'].split("/")[0])
		if score >= score_threshold:
			return True
		else:
			return False
		
	filtered_mash_hits = list(filter(
		lambda x: score_above_threshold(x, mash_hits_score_threshold),
		mash_hits))
	
	# parse plasmid mash
	pathToMashPlasmidScreenTSV = outputDir + "/qcResult/" + ID + "/" + "mashscreen.plasmid.tsv"
	mash_plasmid_hits = result_parsers.parse_mash_result(pathToMashPlasmidScreenTSV)
	# 'shared_hashes' field is string in format '935/1000'
	# Threshold is 100 below highest numerator (ie. '935/100' -> 835)
	mash_plasmid_hits_score_threshold = int(mash_plasmid_hits[0]['shared_hashes'].split("/")[0]) - 100
	filtered_mash_plasmid_hits = list(filter(
		lambda x: score_above_threshold(x, mash_plasmid_hits_score_threshold),
		mash_plasmid_hits))
	
	# parse fastqc
	pathToFastQCR1 = outputDir + "/qcResult/" + ID + "/" + R1[R1.find(os.path.basename(R1)):R1.find(".")] + "_fastqc/summary.txt"
	pathToFastQCR2 = outputDir + "/qcResult/" + ID + "/" + R2[R2.find(os.path.basename(R2)):R2.find(".")] + "_fastqc/summary.txt"
	fastqcR1 = result_parsers.parse_fastqc_result(pathToFastQCR1)
	fastqcR2 = result_parsers.parse_fastqc_result(pathToFastQCR2)
	fastqc = {}
	fastqc["R1"]=fastqcR1
	fastqc["R2"]=fastqcR2
	 
	#parse kraken2 result
	#pathToKrakenResult = outputDir + "/qcResult/" + ID + "/kraken2.genome.report"
	#kraken_genomes = result_parsers.parse_kraken_result(pathToKrakenResult)

	

	multipleSpecies = False
	correctSpecies = False
	containsPlasmid = False
	acceptableCoverage = False
	acceptableFastQCForward = False
	acceptableFastQCReverse = False

	#all the qC result are parsed now, lets do some QC logic
	#look at mash results first
	if (len(filtered_mash_hits) > 1):
		multipleSpecies = True
	
	for mash_hit in filtered_mash_hits:
		species = mash_hit['query_comment']
		if (species.find(expectedSpecies) > -1):
			correctSpecies=True
	
	if (len(filtered_mash_plasmid_hits) > 0):
		containsPlasmid = True
	#look at kraken results > Removed, simplified to using only mash. Validation with Kraken2 is no longer required as we are no longer doing contamination filtering.
	# TODO: filter kraken_genomes based on contamination thresholds
	# if (len(kraken_genomes) > 1):
	#	 output.append("!!!Kraken2 predicted multiple species, possible contamination?")
	#	 notes.append("Kraken2: multiple species, possible contamination.")
	#	 #multipleSpecies = True
	# elif (len(kraken_genomes) == 1):
	#	 multipleSpecies = False
		
	#for kraken_genome in  kraken_genomes:
	#	if (kraken_genome['taxon_name'] == expectedSpecies):
	#		correctSpecies = True

	#look at fastqc results
  
	if (fastqc["R1"]["basic_statistics"] == "PASS" and fastqc["R1"]["per_base_sequence_quality"] == "PASS" and fastqc["R1"]["sequence_length_distribution"] == "PASS" ):
		acceptableFastQCForward = True
	if (fastqc["R2"]["basic_statistics"] == "PASS" and fastqc["R2"]["per_base_sequence_quality"] == "PASS" and fastqc["R2"]["sequence_length_distribution"] == "PASS" ):
		acceptableFastQCReverse = True
	
	#download a reference genome
	print("Downloading reference genomes")
	referenceGenomes = []
	if (not multipleSpecies):
		for mash_hit in filtered_mash_hits: #for all the mash hits, aka reference genomes
			qID = mash_hit['query_id'] #hit genome within mash results
			species_name_start = int(mash_hit['query_comment'].index(".")) + 3 #find the start of species name within query_comment column
			species_name_stop = int (mash_hit['query_comment'].index(",")) #find the end of the species name within query_comment column
			if (mash_hit['query_comment'].find("phiX") > -1):
				species = "PhiX" #phix
			else:
				species = str(mash_hit['query_comment'])[species_name_start: species_name_stop] #assign proper species name for reference genome file name
				# find gcf accession
				# TODO: document this or clean it up to be more readable
				gcf = (qID[:qID.find("_",5)]).replace("_","") #find the full gcf accession for ncbi FTP
				gcf = [gcf[i:i+3] for i in range(0, len(gcf), 3)] #break the gcf accession into k=3

				assembly = qID[:qID.find("_genomic.fna.gz")] #find the assembly name

				#build the urls
				fastaURL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/" + gcf[0] + "/" + gcf[1] + "/" + gcf[2] + "/" + gcf[3] + "/" + assembly + "/" + qID #url to fasta
				assemblyStatURL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/" + gcf[0] + "/" + gcf[1] + "/" + gcf[2] + "/" + gcf[3] + "/" + assembly + "/" + qID + "_assembly_stats.txt" #url to assembly stat
				referencePath = os.path.abspath(outputDir + "/qcResult/" + ID + "/" + species.replace(" ",""))
				referenceGenomes.append(referencePath)

				httpGetFile(url, referencePath + ".gz") #fetch the fasta gz
				httpGetFile(url, referencePath + "_genomeStats.txt") # fetch the genome stat
				with gzip.open(referencePath + ".gz", 'rb') as f:
					gzContent = f.read()
				with open(referencePath, 'wb') as out:
					out.write(gzContent)
				os.remove(referencePath + ".gz")
	else: #throw an error if it contains contaminations
		print("Contaminated Genome assembly...resequencing required")
		raise Exception("contamination and mislabeling...crashing")
		
	#check to make sure we ONLY have ONE reference.
	if (len(referenceGenomes) > 1 ):
		raise Exception ("there are multiple reference genomes")
	elif (len(referenceGenomes) == 0)
		raise Exception ("no reference genoe identified")
		
	
	#now we estimate our coverage using total reads and expected genome size
	
	#find expected genome size
	pathToReferenceStats = referenceGenomes[0] + "_genomeStats.txt"
	with open( pathToReferenceStats, 'r') as referenceStats:
		if (referenceStats.readline().find("all	all	all	all	total-length") > -1): #find the total length stat
			expectedGenomeSize = float(s.split("\t")[5].strip()) 

	#find total base count
	pathToTotalBP = outputDir + "/qcResult/" + ID + "/" + "totalbp"
	with open(pathToTotalBP, 'r') as totalbp_file:
		total_bp = float(totalbp_file.readline())
	
	#calculate coverage
	coverage = total_bp / expectedGenomeSize
			
	if (coverage >= 30):
		acceptableCoverage = True
	
	print(str(multiple))
	print(str(correctSpecies))
	print(str(containsPlasmid))
	print(str(acceptableCoverage))
	print(str(acceptableFastQCForward))
	print(str(acceptableFastQCReverse))

	
	print("Formatting the QC results")

	output.append("\n\n~~~~~~~QC summary~~~~~~~")
	output.append("Estimated genome size: " + str(stats['estimated_genome_size']))
	output.append("Estimated coverage: " + str(stats['estimated_depth_of_coverage']))
	output.append("Expected isolate species: " + expectedSpecies)

	output.append("\nFastQC summary:")
	output.append("\nforward read qc:")
	for key, value in fastqcR1.items():
		output.append(key + ": " + value)
		if (value == "WARN" or value == "FAIL"):
			notes.append("FastQC: Forward read, " + key + " " + value)

	output.append("\nreverse read qc:")
	for key, value in fastqcR2.items():
		output.append(key + ": " + value)
		if (value == "WARN" or value == "FAIL"):
			notes.append("FastQC: Reverse read, " + key + " " + value)

	output.append("\nKraken2 predicted species (>1%): ")
	for kraken_genome in kraken_genomes:
		# TODO filter by species and percentage
		output.append(kraken_genome['taxon_name'])

	output.append("\nmash predicted genomes")
	for mash_hit in filtered_mash_hits:
		output.append(mash_hit['query_comment'])

	output.append("\nmash predicted plasmids")
	for mash_plasmid_hit in mash_plasmid_hits:
		output.append(mash_plasmid_hit['query_comment'])
	
	output.append("\nDetailed kraken genome hits: ")
	for kraken_genome in kraken_genomes:
		output.append(
			str(kraken_genome['fragment_percent']) + '\t' +
			str(kraken_genome['fragment_count_root']) + '\t' +
			str(kraken_genome['fragment_count_taxon']) + '\t' +
			kraken_genome['rank_code'] + '\t' +
			kraken_genome['ncbi_taxon_id'] + '\t' +
			kraken_genome['taxon_name']
		)

	output.append("\nDetailed mash genome hits: ")
	for mash_hit in mash_hits:
		output.append(
			str(mash_hit['identity']) + '\t' +
			mash_hit['shared_hashes'] + '\t' +
			str(mash_hit['median_multiplicity']) + '\t' +
			str(mash_hit['p_value']) + '\t' +
			mash_hit['query_id'] + '\t' +
			mash_hit['query_comment']
		)
	
	output.append("\nDetailed mash plasmid hits: ")
	for mash_plasmid_hit in mash_plasmid_hits:
		output.append(
			str(mash_plasmid_hit['identity']) + '\t' +
			mash_plasmid_hit['shared_hashes'] + '\t' +
			str(mash_plasmid_hit['median_multiplicity']) + '\t' +
			str(mash_plasmid_hit['p_value']) + '\t' +
			mash_plasmid_hit['query_id'] + '\t' +
			mash_plasmid_hit['query_comment']
		)

	#qcsummary
	output.append("\n\nQC Information:")


	'''
	if present:
		output.append("The expected species is predicted by kraken 2")
		#correctSpecies = True
	else:
		output.append("!!!The expected species is NOT predicted by kraken2, contamination? mislabeling?")
		notes.append("Kraken2: Not expected species. Possible contamination or mislabeling")
	
	if (stats['estimated_depth_of_coverage'] < 30):
		output.append("!!!Coverage is lower than 30. Estimated depth: " + str(stats['estimated_depth_of_coverage']))
	
	
	if (len(filtered_mash_hits) > 1):
		output.append("!!!MASH predicted multiple species, possible contamination?")
		multipleSpecies = True
	elif (len(filtered_mash_hits) < 1):
		output.append("!!!MASH had no hits, this is an unknown species")
		notes.append("Mash: Unknown Species")
	'''
	'''
	present=False
	for mash_hit in filtered_mash_hits:
		species = mash_hit['query_comment']
		if (species.find(expectedSpecies) > -1):
			present=True
	
	if present:
		output.append("The expected species is predicted by mash")
		correctSpecies = True
		#notes.append("The expected species is predicted by mash")
	else:
		output.append("!!!The expected species is NOT predicted by mash, poor resolution? contamination? mislabeling?")
		notes.append("Mash: Not expected species. Possible resolution issues, contamination or mislabeling")

	if (len(mash_plasmid_hits) == 0):
		output.append("!!!no plasmids predicted")
		notes.append("Mash: no plasmid predicted")

	#hack: throw exception if this analysis should not proceed due to contamination and mislabelling
	#if (multipleSpecies and not correctSpecies):
	#	out = open(outputDir + "/summary/" + ID +".err", 'a')
	#	out.write('GO RESEQUENCE THIS SAMPLE: It has contamination issues AND mislabelled')
	#	#raise Exception('GO RESEQUENCE THIS SAMPLE: It has contamination issues AND mislabelled')
	
	#identify and download reference genomes, this should only return one reference.
	'''
	
	'''
	correctReference = ""
	if (len(referenceGenomes) > 1):
		for item in referenceGenomes:
			if (item.find(expectedSpecies.replace(" ", "")) > -1): #found the right genome
				correctReference = os.path.basename(item)
	else:
		correctReference = os.path.basename(referenceGenomes[0])
		
	if (correctReference == ""):
		raise Exception("no reference genome...crashing")
	'''
	print("step 2: genome assembly and QC")

	#run the assembly shell script.
	#input parameters: 1 = id, 2= forward, 3 = reverse, 4 = output, 5=tmpdir for shovill, 6=reference genome, 7=buscoDB
	cmd = [scriptDir + "/pipeline_assembly.sh", ID, R1, R2, outputDir, tempDir, referenceGenomes[0], buscodb]
	result = execute(cmd, curDir)
	
	print("Parsing assembly results")
	#get the correct busco and quast result file to parse
	
	'''
	correctAssembly = ""
	buscoPath = "" 
	quastPath = ""
	if ((not multiple and correctSpecies) or (not multiple and not correctSpecies)): #only 1 reference genome
		buscoPath = (outputDir + "/assembly_qc/" + ID + "/" + ID + ".busco" + "/short_summary_" + ID + ".busco.txt")
		quastPath = (outputDir + "/assembly_qc/" + ID + "/" + ID + ".quast" + "/report.tsv")
		correctAssembly = ID
	elif(multiple and correctSpecies): #multiple reference genome, need to find the one we care about
		for item in referenceGenomes:
			if (item.find(expectedSpecies.replace(" ", "")) > -1): #found the right genome
				buscoPath = (outputDir + "/assembly_qc/" + ID + "/" + ID + "." + os.path.basename(item) + ".busco" + "/short_summary_" + ID + "." + os.path.basename(item) + ".busco.txt")
				quastPath = (outputDir + "/assembly_qc/" + ID + "/" + ID + "." + os.path.basename(item) + ".quast/runs_per_reference/" + os.path.basename(item) + "/report.tsv")
				correctAssembly = ID + "." + os.path.basename(item)
	elif(multiple and not correctSpecies):
		raise Exception('GO RESEQUENCE THIS SAMPLE: It has contamination issues AND mislabelled')

	if (buscoPath == "" or quastPath == ""):
		raise Exception('theres no reference genome for this sample for whatever reason...')
	'''
	buscoPath = (outputDir + "/assembly_qc/" + ID + "/" + ID + ".busco" + "/short_summary_" + ID + ".busco.txt")
	quastPath = (outputDir + "/assembly_qc/" + ID + "/" + ID + ".quast" + "/report.tsv")
	#populate the busco and quast result object
	buscoResults = result_parsers.parse_busco_result(buscoPath)
	quastResults = result_parsers.parse_busco_result(quastPath)

if __name__ == "__main__":
	start = time.time()
	print("Starting workflow...")
	main()
	end = time.time()
	print("Finished!\nThe analysis used: " + str(end - start) + " seconds")
