version 1.0

workflow dada2_denoising_iseq {
	input {
		String path_to_fq 
		File path_to_flist
		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"
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
		File path_to_snv
		String no_ref = 'False'
		File reference
		String adjust_mode = "absolute"
		File reference2
		String strain = "3D7"
		String strain2 = "DD2"
		String polyN = "5"
		String min_reads = "0"
		String min_samples = "0"
		String max_snv_dist = "-1"
		String max_indel_dist = "-1"
		String include_failed = "False"
		String exclude_bimeras = "False"
	}
	call ampseq_dada2_iseq_process {
		input:
			path_to_fq = path_to_fq,
			path_to_flist = path_to_flist,
			pattern_fw = pattern_fw,
			pattern_rv = pattern_rv,
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
			overlap_pr2 = overlap_pr2,
			path_to_snv = path_to_snv,
			no_ref = no_ref,
			reference = reference,
			adjust_mode = adjust_mode,
			reference2 = reference2,
			strain = strain,
			strain2 = strain2,
			polyN = polyN,
			min_reads = min_reads,
			min_samples = min_samples,
			max_snv_dist = max_snv_dist,
			max_indel_dist = max_indel_dist,
			include_failed = include_failed,
			exclude_bimeras = exclude_bimeras
	}

	output {
		File rawfilelist_f = ampseq_dada2_iseq_process.rawfilelist
		File missing_files_f = ampseq_dada2_iseq_process.missing_files
		File ASVBimeras_f = ampseq_dada2_iseq_process.ASVBimeras
		File CIGARVariants_Bfilter_f = ampseq_dada2_iseq_process.CIGARVariants_Bfilter
		File ASV_to_CIGAR_f = ampseq_dada2_iseq_process.ASV_to_CIGAR
		File seqtab_f = ampseq_dada2_iseq_process.seqtab
	}
}

task ampseq_dada2_iseq_process {
	input {
		String path_to_fq 
		File path_to_flist
		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"
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
		File path_to_snv
		String no_ref = 'False'
		File reference
		String adjust_mode = "absolute"
		File reference2
		String strain = "3D7"
		String strain2 = "DD2"
		String polyN = "5"
		String min_reads = "0"
		String min_samples = "0"
		String max_snv_dist = "-1"
		String max_indel_dist = "-1"
		String include_failed = "False"
		String exclude_bimeras = "False"
	}

	Map[String, String] in_map = {
		"path_to_fq": "fq_dir",
		"path_to_flist": sub(path_to_flist, "gs://", "/cromwell_root/"),
		"pattern_fw": pattern_fw,
		"pattern_rv": pattern_rv,
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
		"overlap_pr2" : sub(overlap_pr2, "gs://", "/cromwell_root/"),
		"path_to_snv": sub(path_to_snv, "gs://", "/cromwell_root/"),
		"no_ref": no_ref,
		"reference": sub(reference, "gs://", "/cromwell_root/"),
		"adjust_mode": adjust_mode,
		"reference2": sub(reference2, "gs://", "/cromwell_root/"),
		"strain": strain,
		"strain2": strain2,
		"polyN": polyN,
		"min_reads": min_reads,
		"min_samples": min_samples,
		"max_snv_dist": max_snv_dist,
		"max_indel_dist": max_indel_dist,
		"include_failed": include_failed,
		"exclude_bimeras": exclude_bimeras
	}
	File config_json = write_json(in_map)
	command <<<
	set -euxo pipefail
	#set -x
	mkdir fq_dir

	gsutil ls ~{path_to_fq}
	gsutil -m cp -r ~{path_to_fq}* fq_dir/

	python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --mixed_reads --meta --repo --adaptor_removal --dada2
	cat Results/stderr.txt
	cat Results/stdout.txt
	cat Results/AdaptorRem/*stderr.txt
	cat Results/AdaptorRem/adaptorrem_meta.tsv

	find . -type f
	python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --mixed_reads --dada2 #--postproc_dada2 --asv_to_cigar
	cat Results/stderr.txt
	cat Results/stdout.txt
	find . -type f
	>>>
	output {
		File rawfilelist = "Results/Fq_metadata/rawfilelist.tsv"
		File missing_files = "Results/missing_files.tsv" 
		File ASVBimeras = "Results/ASVBimeras.txt"
		File CIGARVariants_Bfilter = "Results/CIGARVariants_Bfilter.out.tsv"
		File ASV_to_CIGAR = "Results/ASV_to_CIGAR/ASV_to_CIGAR.out.txt"
		File seqtab = "Results/seqtab_mixed.tsv"
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
