"""\
Amplicon decontamination scripts
"""

import pandas as pd
import os, fnmatch
import subprocess
import sys
import shutil
import glob
import gzip

from Bio import SeqIO
from Bio.Seq import reverse_complement
from Bio.Seq import Seq

#GENERAL USE SCRIPTS SECTION
def flush_dir(parent_dir, *dirnames):
	"""
	Remove all files and subdirectories from a directory, and create a new empty
	directory with the same name. Multiple subdirectiories may be provided.

	Args:
	parent_dir (str): The path of the parent directory.
	dirname (str): The name of the directory to be flushed.

	Returns:
	None
	"""

	dirpath = os.path.join(parent_dir, *dirnames)
	shutil.rmtree(dirpath, ignore_errors=True)
	os.makedirs(dirpath)
	return ()

def create_meta(path_to_fq, parent_dir, dirname, filename, pattern_fw, pattern_rv):
	"""
	Creates a metadata file with the list of files to process and their paths.

	Args:
	- path_to_fq: string, path to the directory containing the input fastq files
	- parent_dir: string, path to the parent directory where the output files will be stored
	- dirname: string, name of the subdirectory where the output files will be stored
	- filename: string, name of the output metadata file
	- pattern_fw: string, pattern to match forward reads
	- pattern_rv: string, pattern to match reverse reads

	Returns:
	- None

	Example usage:

	create_meta('/path/to/input/fastq/', '/path/to/output/', 'subdir', 'metadata.tsv', '*_R1.fastq', '*_R2.fastq')

	This will search for all files in /path/to/input/fastq/ that end with '_R1.fastq' and '_R2.fastq',
	create a metadata file named 'metadata.tsv', and store it in a subdirectory named 'subdir' within 
	/path/to/output/.
	"""

	filelist = os.listdir(path_to_fq)

	dirpath = os.path.join(parent_dir, dirname) 
	outfile = os.path.join(dirpath, filename)
	print(f"Meta file will be generated at {outfile}")

	meta_df = []
	for file_fw in glob.glob(os.path.join(path_to_fq, pattern_fw)):
		sampleid = os.path.basename(file_fw).split(pattern_fw[1:])[0]
		file_rv = os.path.join(path_to_fq, sampleid + pattern_rv[1:])
		if os.path.isfile(file_rv):
			meta_df.append((sampleid, file_fw, file_rv))

	with open(outfile, 'w') as f:
		for row in meta_df:
			f.write('\t'.join(row) + '\n')

	print(f"Meta file generated at {outfile}")
	return()

#MERGE FUNCTIONS SECTION
def mergereads(sampleid, fileF, fileR, res_dir, subdir, read_maxlength=200, pairread_minlength=100, merge_minlength=100):
	"""
	This function uses bbmerge.sh to merge paired-end reads from two fastq files
	(fileF and fileR) into a single fastq file. It also generates two other fastq
	files of unmerged reads. The output files are saved in the specified res_dir
	and subdir directory paths. The function also creates a metadata file 
	(merge_meta.tsv) containing the sample ID, output file name, and standard 
	output and error logs. If either the forward or reverse fastq files are not 
	found, the function exits with an error message. This functions is optimized 
	for reads that are at the least 200bp long and amplicons 100 bp long or longer.
	Merging shorter reads will require chaning this parameters in the config file.
	
	Args:
	sampleid: a string representing the sample identifier.
	fileF: a string representing the file path of the forward reads.
	fileR: a string representing the file path of the reverse reads.
	res_dir: a string representing the directory path of the results.
	subdir: a string representing the subdirectory path of the results.
	read_maxlength: an integer representing the maximum length of the read. Read after this will be trimmed.
	pairread_minlength: an integer representing the minimum length of the mated reads.
	merge_minlength: an integer representing the minimum merge length. Merge shorter than this will be discarded.

	Returns: None

	Example usage:

	mergereads("sample1", "/path/to/forward.fastq", "/path/to/reverse.fastq", "/path/to/results", "subdirectory")
	"""

	if os.path.isfile(fileF) and os.path.isfile(fileR):	
		file_nameout = os.path.join(res_dir, subdir, f"{sampleid}_stdout.txt")
		file_nameerr = os.path.join(res_dir, subdir, f"{sampleid}_stderr.txt")
		output_file_path = os.path.join(res_dir, subdir, f"{sampleid}_merged.fastq")
		output_unmerged_f_path = os.path.join(res_dir, subdir, f"{sampleid}_unmergedf.fastq")
		output_unmerged_r_path = os.path.join(res_dir, subdir, f"{sampleid}_unmergedr.fastq")
		meta_file_path = os.path.join(res_dir, subdir, f"merge_meta.tsv")

		sys.stdout = open(file_nameout, "w")
		sys.stderr = open(file_nameerr, "w")

		with open(meta_file_path, "a") as meta_file:
			meta_file.write(f"{sampleid}\t{output_file_path}\t{file_nameout}\t{file_nameerr}\n")
		
		#The meaning of insert in mininsert is different to the meaning of insert in the cont_report.
		#mininsert equals the final size of the merged read. That is, fbarcode + insert + rbarcode in cont_report. 	
		cmd = ['bbmerge.sh', 
			f'forcetrimright={read_maxlength}', 
			f'minlength={pairread_minlength}', 
			f'mininsert={merge_minlength}', 
			f'in1={fileF}', 
			f'in2={fileR}', 
			f'out={output_file_path}',
			f'outu1={output_unmerged_f_path}', 
			f'outu2={output_unmerged_r_path}']
		proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
		proc.wait()

		sys.stdout = open(os.path.join(res_dir, 'stdout.txt'), 'a')
		sys.stderr = open(os.path.join(res_dir, 'stderr.txt'), 'a')

	else:
		sys.exit('BBmerge halted : one or both of the fastq files not found! Exiting..')
	return()

def extract_bbmergefields(sampleid, mergefile, bbreportfile, path_to_flist, res_dir, rep_dir, subdir):
	"""
	Extracts relevant data from a bbmerge report file and saves it to a tab-separated file.

	Args:
	sampleid: the ID of the sample being processed
	mergefile: the path to the file with the merged reads
	bbreportfile: the path to the bbmerge report file
	path_to_flist: the path to a csv file with the sample_id,Forward,Reverse, where Forward and Reverse are columns with the barcodes for the sample
	res_dir: the path to the main results directory
	rep_dir: the path to the reports directory within the results directory
	subdir: the name of the subdirectory within the results directory where output files should be written

	Returns:
	None

	Example Usage:

	extract_bbmergefields("Sample1", "/path/to/bbmerge.fastq", "/path/to/bbmerge_report.txt", "path/to/barcodes_match.csv", "/path/to/results", "/path/to/reports", "bbmerge")
	"""

	if os.path.isfile(bbreportfile) and os.path.isfile(mergefile):				
		sys.stdout = open(os.path.join(rep_dir, 'stdout.txt'), 'a')
		sys.stderr = open(os.path.join(rep_dir, 'stderr.txt'), 'a')
		bbmergedata = {}
		bbmergedata['sampleid'] = sampleid

		with open(bbreportfile, 'r') as f, open(os.path.join(rep_dir, subdir, "bbmergefields.tsv"), 'a') as o:
			for line in f:
				fields = line.split()
				if not fields:	#Skip empty lines
					continue	
				if fields[0] == 'Pairs:':
					bbmergedata['Pairs'] = fields[1]
				elif fields[0] == 'Joined:':
					bbmergedata['Joined'] = fields[1]
					bbmergedata['JoinedP'] = fields[2].strip('%')
				elif fields[0] == 'Ambiguous:':
					bbmergedata['Ambiguous'] = fields[1]
					bbmergedata['AmbiguousP'] = fields[2].strip('%') 
				elif fields[0] == 'No':
					bbmergedata['No_Solution'] = fields[2]
					bbmergedata['No_SolutionP'] = fields[3].strip('%') 
				elif fields[0] == 'Too':
					bbmergedata['Too_Short'] = fields[2]
					bbmergedata['Too_ShortP'] = fields[3].strip('%') 
				elif fields[0] == 'Avg':
					bbmergedata['Avg_Insert'] = fields[2]
				elif fields[0] == 'Standard':
					bbmergedata['Standard_Deviation'] = fields[2]
				elif fields[0] == 'Mode:':
					bbmergedata['Mode'] = fields[1]
				elif fields[0] == 'Insert':
					bbmergedata['Insert_range_low'] = fields[2]
					bbmergedata['Insert_range_up'] = fields[4]
				elif fields[0] == '90th':
					bbmergedata['90th_pc'] = fields[2]
				elif fields[0] == '75th':
					bbmergedata['75th'] = fields[2]
				elif fields[0] == '50th':
					bbmergedata['50th_pc'] = fields[2]
				elif fields[0] == '25th':
					bbmergedata['25th_pc'] = fields[2]
				elif fields[0] == '10th':
					bbmergedata['10th_pc'] = fields[2]

			o.write('\t'.join(bbmergedata.values()) + '\n')

		cmd = ['Rscript', os.path.join('/Code/runBBMergecontamination.R'),
		'-p', f'{mergefile}',
		'-d', os.path.join(rep_dir, subdir),
		'-b', path_to_flist]

		print(cmd)
		proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
		proc.wait()

		sys.stdout = open(os.path.join(res_dir, 'stdout.txt'), 'a')
		sys.stderr = open(os.path.join(res_dir, 'stderr.txt'), 'a')
	else:
		sys.exit('Extract bbmerge report halted : bbmerge report file not found! Exiting..')
	return()

#ADAPTOR AND PRIMER REMOVAL SECTION						
def adaptor_rem(sampleid, fileF, fileR, res_dir, subdir, qvalue = 5, length = 20): 
	"""
	Runs the Trim Galore tool to remove adaptors and trim low-quality reads from paired-end fastq files.

	Args:
	sampleid (str): The base name for the output files.
	fileF (str): The path to the forward-read fastq file.
	fileR (str): The path to the reverse-read fastq file.
	res_dir (str): The path to the directory where the output files will be saved.
	subdir (str): The name of the subdirectory within the results directory where output files should be written
	qvalue (int, optional): The minimum quality score for trimming. Defaults to 5.
	length (int, optional): The minimum length of the reads to keep after trimming. Defaults to 20.

	Returns:
	None
	"""

	if os.path.isfile(fileF) and os.path.isfile(fileR):
		file_nameout = os.path.join(res_dir, subdir, f"{sampleid}_stdout.txt")
		file_nameerr = os.path.join(res_dir, subdir, f"{sampleid}_stderr.txt")
		sys.stdout = open(file_nameout, "w")
		sys.stderr = open(file_nameerr, "w")

		output_dir = os.path.join(res_dir, subdir)
		cmd = ['trim_galore', '--paired', '--gzip', '--quality', f'{qvalue}', '--length', 
		f'{length}', '--output_dir', f'{output_dir}', '--basename', f'{sampleid}', f'{fileF}', 
		f'{fileR}']
		proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
		#proc = subprocess.Popen(cmd)#, stdout=sys.stdout, stderr=sys.stderr)
		proc.wait()

		sys.stdout = open(os.path.join(res_dir, 'stdout.txt'), 'a')
		sys.stderr = open(os.path.join(res_dir, 'stderr.txt'), 'a')
	else:
		sys.exit('Adaptor Removal halted : one or both of the fastq files not found! Exiting..')
	return()

def trim_primer(sampleid, fileF, fileR, res_dir, subdir, pr1, pr2, prefix, keep_untrimmed=False):
	"""
	Trim primers from paired-end fastq files using cutadapt.

	Args:
	sampleid (str): Sample identifier.
	fileF (str): Path to input forward fastq file.
	fileR (str): Path to input reverse fastq file.
	res_dir (str): Path to output directory.
	subdir (str): The name of the subdirectory within the results directory where output files should be written
	pr1 (str): Path to primer sequence file for forward read.
	pr2 (str): Path to primer sequence file for reverse read.
	prefix (str): Prefix to use for output filenames.
	keep_untrimmed (bool, optional): If True, keep untrimmed reads in separate files. Default is False.

	Returns:
	None
	"""

	if os.path.isfile(fileF) and os.path.isfile(fileR):
		file_nameout = os.path.join(res_dir, subdir, f"{sampleid}_stdout.txt")
		file_nameerr = os.path.join(res_dir, subdir, f"{sampleid}_stderr.txt")
		sys.stdout = open(file_nameout, "w")
		sys.stderr = open(file_nameerr, "w")

		cmd = ['cutadapt', '-g', f'file:{pr1}', '-G', f'file:{pr2}',
			'-o', os.path.join(res_dir, subdir, f'{sampleid}_{prefix}_1.fq.gz'),
			'-p', os.path.join(res_dir, subdir, f'{sampleid}_{prefix}_2.fq.gz'),
			'--pair-adapters', '--action=trim']

		if keep_untrimmed:
			cmd.extend(['--untrimmed-output', os.path.join(res_dir, subdir, f'{sampleid}_temp_1.fq.gz'),
	    			'--untrimmed-paired-output', os.path.join(res_dir, subdir, f'{sampleid}_temp_2.fq.gz')])
		else:
			cmd.append('--discard-untrimmed')

		cmd.extend([fileF, fileR])
		print(cmd)
		proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
		proc.wait()

		sys.stdout = open(os.path.join(res_dir, 'stdout.txt'), 'a')
		sys.stderr = open(os.path.join(res_dir, 'stderr.txt'), 'a')
	else:
		sys.exit('Pre-process halted : one or both of the fastq files not found! Exiting..')
	return()

#RUN DADA2 SECTION
def run_dada2(path_to_DADA2, path_to_meta, path_to_fq, path_to_flist, Class, maxEE, trimRight, minLen, truncQ, matchIDs, max_consist, omegaA, justConcatenate, maxMismatch, saveRdata, res_dir, subdir):
	"""
	Runs the DADA2 pipeline on the input files using the specified parameters.

	Args:
	- path_to_DADA2 (str): the path to the DADA2 installation directory.
	- path_to_meta (str): the path to the metadata file containing sample information.
	- path_to_fq (str): the path to the raw fastq.gz files.
	- path_to_flist (str): the path to a csv file with the sample_id,Forward,Reverse, where Forward and Reverse are columns with the barcodes for the sample
	- res_dir (str): the path to the directory where results will be saved.
	- subdir (str): the name of the subdirectory where the output files will be saved.
	- Class (str): the name of the column in the metadata file that contains the sample class information.
	- maxEE (float): the maximum expected error rate.
	- trimRight (int): the number of bases to trim from the 3' end of the reads.
	- minLen (int): the minimum length of reads to retain after trimming.
	- truncQ (int): the quality threshold for truncating reads.
	- matchIDs (bool): boolean to request DADA2 to match ids on fastqs to make sure reads on forward and reverse end are in same order.
	- max_consist (int): the maximum number of mismatches allowed in the overlap region for merging paired-end reads.
	- omegaA (float): the alpha parameter for the consensus quality score.
	- justConcatenate (int): whether to just concatenate the forward and reverse reads without merging them.
	- maxMismatch (int): the maximum number of mismatches allowed during merging.
	- saveRdata (str): whether to save the intermediate R data files.

	Returns:
	- None
	"""

	if os.path.isfile(path_to_meta):
		file_nameout = os.path.join(res_dir, subdir, "stdout.txt")
		file_nameerr = os.path.join(res_dir, subdir, "stderr.txt")
		sys.stdout = open(file_nameout, "w")
		sys.stderr = open(file_nameerr, "w")

		if any(subdir == i for i in ['DADA2', 'DADA2_OP', 'DADA2_NOP']):
			program = 'runDADA2.R'
		else:
			program = 'runDADA2contamination.R' 

		bimera = '--bimera'
		cmd = ['Rscript', os.path.join("/", path_to_DADA2, program),
		'-p', f'{path_to_meta}',
		'-r', f'{path_to_fq}',
		'-d', os.path.join(res_dir, subdir),
		'-o', os.path.join(res_dir, subdir, 'seqtab.tsv'),
		'-c', f'{Class}',
		'-ee', f'{maxEE}',
		'-tR', f'{trimRight}',
		'-mL', f'{minLen}',
		'-tQ', f'{truncQ}',
		'-id', f'{matchIDs}',
		'-mC', f'{max_consist}',
		'-wA', f'{omegaA}',
		'-jC', f'{justConcatenate}',
		'-mM', f'{maxMismatch}',
		'-s', f'{saveRdata}',
		'-b', f'{path_to_flist}',
		f'{bimera}']
		print(cmd)
		proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
		proc.wait()

		sys.stdout = open(os.path.join(res_dir, 'stdout.txt'), 'a')
		sys.stderr = open(os.path.join(res_dir, 'stderr.txt'), 'a')
	else:
		sys.exit('DADA2 halted : No path to meta file provided! Exiting..')
	return()

def merge_seqtab(path_op, path_nop):
	"""
	Merges overlapping and non-overlapping dada2 tables into a single table.

	Parameters:
	path_op (str): File path to the overlapping dada2 table in CSV format.
	path_nop (str): File path to the non-overlapping dada2 table in CSV format.

	Returns:
	seqtab (pd.DataFrame): Merged dada2 table containing both overlapping and non-overlapping data.

	"""

	if os.path.isfile(path_op) and os.path.isfile(path_nop):
		seqtab_op = pd.read_csv(path_op, sep = "\t")
		seqtab_nop = pd.read_csv(path_nop, sep = "\t")
		seqtab = pd.concat([seqtab_op,seqtab_nop],axis=1)
	else:
		sys.exit('Overlapping and/or non-overlapping dada2 tables not found! Exiting..')

	return(seqtab)

