
#find ./ -name "sketch" -exec rm -r "{}" \;
#find ./ -name "distance" -exec rm -r "{}" \;

import sys, os, shutil, glob, gzip
os.chdir(os.path.split(os.path.realpath(__file__))[0])

suffix = " > /dev/null 2>&1"
dir = "__results__"

K = [21, 31]
FILTER = ["", "-filter"]
NB_READS = ["0", "100"]
NB_KMERS = ["100", "1000"]
NB_CORES = [1, 0]

def create_command(outputPrefix, k, filter, nb_reads, nb_kmers, nb_cores):
	outputDir = "k" + str(k) + "_" + filter.replace("-", "") + "_" + nb_reads + "-" + nb_kmers + "_n" + str(nb_cores)
	command = "python ../../scripts/simkaMin/simkaMin_pipeline.py "
	command += " -in " + " ../../example/simka_input.txt "
	command += " -out " + outputPrefix + "/" + outputDir
	command += " -nb-cores " + str(nb_cores)
	command += " -max-memory 100 "
	command += " -kmer-size " + str(k)
	command += " -nb-kmers " + str(nb_kmers)
	command += " -bin ../../build/bin/simkaMin "
	command += " -max-reads " + str(nb_reads)
	command += " " + filter + " "
	return command, outputDir

def create_truth():
	for k in K:
		for filter in FILTER:
			for nb_reads in NB_READS:
				for nb_kmers in NB_KMERS:
					for nb_cores in NB_CORES:
						command, outputDir = create_command("truth_simkaMin", k, filter, nb_reads, nb_kmers, nb_cores)
						print command
						ret = os.system(command)
						if ret != 0: exit(1)

#create_truth()
#exit(1)

def clear():
	#if os.path.exists("temp_output"):
	#	shutil.rmtree("temp_output")
	if os.path.exists("__results__"):
		shutil.rmtree("__results__")
	os.mkdir(dir)

def decompress_simka_results(dir):
	result_filenames = glob.glob(os.path.join(dir, '*.csv.gz'))
	for filename_gz in result_filenames:
		#filename_gz = result_dir + "/" + filename
		with gzip.open(filename_gz, 'rb') as f:
			outFile = open(filename_gz[:-3], "w")
			outFile.write(f.read())
			outFile.close()
			os.remove(filename_gz)

def __test_matrices(result_dir, truth_dir):

	ok = True

	decompress_simka_results(result_dir)
	result_filenames = glob.glob(os.path.join(result_dir, '*.csv'))
	if len(result_filenames) == 0:
		print("Error: no results")
		exit(1)

	decompress_simka_results(truth_dir)
	truth_filenames = glob.glob(os.path.join(truth_dir, '*.csv'))

	#if simka_vs_truth:
	#	truth_filenames = glob.glob(os.path.join(truth_dir, '*.csv'))
	#else: #simka vs simka
	#	#if result_dir+"/mat_abundance_jaccard.csv" in truth_filenames: #comparing simka results vs simka results
	#	#truth_filenames.remove(result_dir+"/mat_abundance_jaccard.csv") #This distance is computed from Bray Curtis distance
	#
	#	truth_filenames = glob.glob(os.path.join(truth_dir, '*.csv'))

	truth_filenames.sort()
	result_filenames.sort()

	for result_filename in result_filenames:
		distanceName = os.path.split(result_filename)[1]
		for truth_filename in truth_filenames:
			distanceName2 = os.path.split(truth_filename)[1]
			if distanceName != distanceName2: continue

			res_file = open(result_filename, "r")
			truth_file = open(truth_filename, "r")

			#print res_file, truth_file
			res_str = res_file.read()
			truth_str = truth_file.read()

			res_file.close()
			truth_file.close()

			if(res_str != truth_str):
				print("\t- TEST ERROR:    " + distanceName)
				ok = False

	return ok


def test_dists(dir):
	if(__test_matrices("__results__/" + dir, "truth_simkaMin/" + dir)):
		print("\tOK")
	else:
		print("\tFAILED")
		exit(1)


#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
clear()

def test():
	for k in K:
		for filter in FILTER:
			for nb_reads in NB_READS:
				for nb_kmers in NB_KMERS:
					for nb_cores in NB_CORES:
						command, outputDir = create_command(dir, k, filter, nb_reads, nb_kmers, nb_cores)
						print command
						ret = os.system(command + suffix)
						if ret != 0: exit(1)
						test_dists(outputDir)
						clear()



#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
def test_append():
	print("Test append command")
	clear()

	out_dir = "./test_append"
	if os.path.exists(out_dir):
		shutil.rmtree(out_dir)
	os.mkdir(out_dir)

	merged_sketch_filename = os.path.join(out_dir, "merged_sketch.bin")
	filename = "../../example/simka_input.txt"
	for line in open(filename):
		line = line.strip()
		if len(line) == 0: continue

		filename_temp = os.path.join("../../example/test_simkaMin_input_temp.txt")
		f = open(filename_temp, "w")
		f.write(line)
		f.close()

		sketch_filename = os.path.join(out_dir, "sketch.bin")
		command = "../../build/bin/simkaMin sketch -in " + filename_temp + " -out " + sketch_filename + " -nb-kmers 100 -kmer-size 21"
		ret = os.system(command + suffix)
		if ret != 0: exit(1)


		if os.path.exists(merged_sketch_filename):
			command = "../../build/bin/simkaMin append -in1 " + merged_sketch_filename + " -in2 " + sketch_filename
			ret = os.system(command + suffix)
			if ret != 0: exit(1)
			os.remove(sketch_filename)
		else:
			shutil.move(sketch_filename, merged_sketch_filename)

		os.remove(filename_temp)

	command = "../../build/bin/simkaMin distance -in1 " +  merged_sketch_filename + " -in2 " + merged_sketch_filename + " -out " + dir
	ret = os.system(command + suffix)
	if ret != 0: exit(1)

	command = "../../build/bin/simkaMin export -in " + dir + " -in1 " +  merged_sketch_filename + " -in2 " + merged_sketch_filename + " -out " + dir
	ret = os.system(command + suffix)
	if ret != 0: exit(1)


	if(__test_matrices(dir, "truth_simkaMin/k21__0-100")):
		print("\tOK")
	else:
		print("\tFAILED")
		exit(1)


#test()
test_append()



if os.path.exists("__results__"):
	shutil.rmtree("__results__")
if os.path.exists("test_append"):
	shutil.rmtree("test_append")