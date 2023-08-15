version 1.0

workflow amplicon_decontamination_detect {
	input {
		String path_to_fq 
		File path_to_flist
		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"
		Int joined_threshold = 1000
		Float contamination_threshold = 0.5
		File pr1 
		File pr2 
		String Class = "parasite"
		String maxEE = "5,5"
		String trimRight = "0,0"
		Int minLen = 30
		String truncQ = "5,5"
		String matchIDs = "0"
		Int max_consist = 10
		Float omegaA = 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
		String saveRdata = ""
		Int justConcatenate = 0
		Int maxMismatch = 0
		String path_to_DADA2 = '/Code'
		File overlap_pr1
		File overlap_pr2
	}
	call ampseq_dada2merge_iseq_process {
		input:
			path_to_fq = path_to_fq,
			path_to_flist = path_to_flist,
			pattern_fw = pattern_fw,
			pattern_rv = pattern_rv,
			joined_threshold = joined_threshold,
			contamination_threshold = contamination_threshold,
			pr1 = pr1,
			pr2 = pr2,
			Class = Class,
			maxEE = maxEE,
			trimRight = trimRight,
			minLen = minLen,
			truncQ = truncQ,
			matchIDs = matchIDs,
			max_consist = max_consist,
			omegaA = omegaA,
			saveRdata = saveRdata,
			justConcatenate = justConcatenate,
			maxMismatch = maxMismatch,
			path_to_DADA2 = path_to_DADA2,
			overlap_pr1 = overlap_pr1,
			overlap_pr2 = overlap_pr2
	}

	output {
		File rawfilelist_f = ampseq_dada2merge_iseq_process.rawfilelist
		File missing_files_f = ampseq_dada2merge_iseq_process.missing_files
		File Barcode_report_abs_f = ampseq_dada2merge_iseq_process.Barcode_report_abs
		File Barcode_report_per_f = ampseq_dada2merge_iseq_process.Barcode_report_per
		File Insert_size_f = ampseq_dada2merge_iseq_process.Insert_size
		File Match_report_abs_f = ampseq_dada2merge_iseq_process.Match_report_abs
		File Match_report_per_f = ampseq_dada2merge_iseq_process.Match_report_per
		File barcodes_report_dada2_f = ampseq_dada2merge_iseq_process.barcodes_report_dada2
		File hamming_distances_forward_f = ampseq_dada2merge_iseq_process.hamming_distances_forward
		File hamming_distances_reverse_f = ampseq_dada2merge_iseq_process.hamming_distances_reverse	
		File seqtab_f = ampseq_dada2merge_iseq_process.seqtab
		File sequences_barplot_f = ampseq_dada2merge_iseq_process.sequences_barplot
		File stacked_barplot_per_f = ampseq_dada2merge_iseq_process.stacked_barplot_per
		File stacked_barplot_f = ampseq_dada2merge_iseq_process.stacked_barplot
	}
}

task ampseq_dada2merge_iseq_process {
	input {
		String path_to_fq 
		File path_to_flist
		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"
		Int joined_threshold = 1000
		Float contamination_threshold = 0.5
		File pr1 
		File pr2
		String Class = "parasite"
		String maxEE = "5,5"
		String trimRight = "0,0"
		Int minLen = 30
		String truncQ = "5,5"
		String matchIDs = "0"
		Int max_consist = 10
		Float omegaA = 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
		String saveRdata = ""
		Int justConcatenate = 0
		Int maxMismatch = 0
		String path_to_DADA2 = '/Code'
		File overlap_pr1
		File overlap_pr2
	}

	Map[String, String] in_map = {
		"path_to_fq": "fq_dir",
		"path_to_flist": sub(path_to_flist, "gs://", "/cromwell_root/"),
		"pattern_fw": pattern_fw,
		"pattern_rv": pattern_rv,
		"joined_threshold": joined_threshold,
		"contamination_threshold": contamination_threshold,
		"pr1": sub(pr1, "gs://", "/cromwell_root/"),
		"pr2": sub(pr2, "gs://", "/cromwell_root/"),
		"Class": Class,
		"maxEE": maxEE,
		"trimRight": trimRight,
		"minLen": minLen,
		"truncQ": truncQ,
		"matchIDs": matchIDs,
		"max_consist": max_consist,
		"omegaA": omegaA,
		"saveRdata": saveRdata,
		"justConcatenate": justConcatenate,
		"maxMismatch": maxMismatch,
		"path_to_DADA2": path_to_DADA2,
		"overlap_pr1" : sub(overlap_pr1, "gs://", "/cromwell_root/"),
		"overlap_pr2" : sub(overlap_pr2, "gs://", "/cromwell_root/")
	}
	File config_json = write_json(in_map)
	command <<<
	set -euxo pipefail
	#set -x
	mkdir fq_dir
	
	gsutil ls ~{path_to_fq}
	gsutil -m cp -r ~{path_to_fq}* fq_dir/

	python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --mixed_reads --meta --repo --adaptor_removal --primer_removal --dada2_contamination

	Rscript /Code/Contamination.R Report/DADA2_Contamination/ Report/ ~{path_to_flist} ~{joined_threshold} ~{contamination_threshold}
	find . -type f
	>>>
	output {
		File rawfilelist = "Results/Fq_metadata/rawfilelist.tsv"
		File missing_files = "Results/missing_files.tsv" 
		File Barcode_report_abs = "Report/Barcode_report_abs.svg"
		File Barcode_report_per = "Report/Barcode_report_per.svg"
		File Insert_size = "Report/Insert_size.png"
		File Match_report_abs = "Report/Match_report_abs.svg"
		File Match_report_per = "Report/Match_report_per.svg"
		File barcodes_report_dada2 = "Report/barcodes_report_dada2.tsv"
		File hamming_distances_forward = "Report/hamming_forward.tsv"
		File hamming_distances_reverse = "Report/hamming_reverse.tsv"
		File seqtab = "Results/DADA2_Contamination/seqtab.tsv"
		File sequences_barplot = "Results/DADA2_Contamination/sequences_barplot.svg"
		File stacked_barplot_per = "Results/DADA2_Contamination/stacked_barplot_per.svg"
		File stacked_barplot = "Results/DADA2_Contamination/stacked_barplot.svg"

	}
	runtime {
		cpu: 1
		memory: "10 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: 'jorgeamaya/dada2_terra_iseq:v1'
	}
}
