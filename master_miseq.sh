#!/bin/bash

###############################################################################
#SCRIPT NAME: Amplicon Terra Pipeline master script    	 		      #
#DESCRIPTION: Analyse Amplicon Sequencing Data for Contamination	      #
#ARGS: No arguments						              #
#AUTHORS: J.E. Amaya Romero, Angela Early, Jason Mohabir, Phillip Schwabl     #	
#CONTACT: jamayaro@broadinstitute.org	       				      #
###############################################################################

:<<'CONFIG.JSON'
#Instructions for the CONFIG.JSON FILE

#path_to_fq: Path to input fastq files
#path_to_flist: Path to list with all the expected samples.
#pattern_fw: Pattern for forward reads in fastqs names
#pattern_rv: Pattern for reverse reads in fastqs names
#read_maxlength: Crop reads at this length. Use to avoid merging bad quality reads.
#pairread_minlength: Minimum paired read length. Use to remove unusual sequence pairs. Amplicons must be longer than this.
#merge_minlength: Minimum merge length. Use to remove unusual sequence pairs. Amplicons must be longer than this. Must be equal or shorter than pairread_minlength.
#barcodes_file: Fasta file with barcodes
#pr1: Path to forward primers FASTA file
#pr2: Path to reverse primers FASTA file
#Class: Specify Analysis class. Accepts one of two: parasite/vector
#maxEE: Maximum Expected errors (dada2 filtering argument)
#trimRight: Hard trim number of bases at 5 end (dada2 filtering argument)
#minLen: Minimum length filter (dada2 filtering argument)
#truncQ: Soft trim bases based on quality (dada2 filtering argument)
#max_consist: Number of cycles for consistency in error model (dada2 argument)
#omegaA: p-value for the partitioning algorithm (dada2 argument)
#saveRdata: Optionally save dada2 part of this run as Rdata object
#justConcatenate: whether reads should be concatenated with N's during merge (dada2 argument)
#maxMismatch: Specify the maximum number of mismatches allowed during merging
#overlap_pr1: Path to forward primers for shorter overlapping targets FASTA file (For mixed_reads run only)
#overlap_pr2: Path to reverse primers for shorter overlapping targets FASTA file (For mixed_reads run only)
#reference: Path to reference target sequences (If --mixed_reads flag is set)
#strain: Main strain for the postprocess of DADA2 ASVs
#strain2: Optional second strain for the postprocess of DADA2 ASVs
#polyN: Mask homopolymer runs length >= polyN (default: 5; disabled < 2)
#min_reads: Minimum total reads to include ASV (default: 0, disabled)
#min_samples: Minimum samples to include ASV (default: 0, disabled)
#max_snv_dist: Maximum SNV distance to include ASV (default: -1, disabled)
#max_indel_dist: Maximum indel distance to include ASV (default: -1, disabled)
#include_failed: INCLUDE ASVs that failed post-DADA2 filters (default: False)
#exclude_bimeras: EXCLUDE ASVs that DADA2 flagged as bimeras (default: False)
#amp_mask: Amplicon low complexity mask info (default: None, disabled)
#verbose: Increase verbosity
CONFIG.JSON

#EXEC

#conda activate ampseq_env

python Code/Amplicon_TerraPipeline.py --config config_MiSeq.json --overlap_reads \
--meta \
--repo \
--adaptor_removal \
--primer_removal \
--dada2_contamination #\
#--dada2 \
#--postproc_dada2 \
#--asv_to_cigar

#if [ -d "$PWD/Report/Merge/" ]; then
	#These plots are optimized for datasets up to 100 specimens. This is by design as a standard illumina plate has 96 wells. Larger dataset will still produce the plots, but some column names will be missing. Names longer than 60 characters will extend outside the plotting area.
#	Rscript Code/BBMerge.R "$PWD/Report/Merge/" "$PWD/Report/"
#	Rscript Code/Contamination.R "$PWD/Report/Merge/" "$PWD/Report/" "$PWD/Data/barcodes_matches.csv" '1000' '0.5'
#fi

#if [ -d "$PWD/Report/DADA2_Contamination/" ]; then
#	Rscript Code/Contamination.R "$PWD/Report/DADA2_Contamination/" "$PWD/Report/" "$PWD/Data/barcodes_matches.csv" '1000' '0.5'
#fi
