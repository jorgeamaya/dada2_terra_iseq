#!/bin/bash

cut -f1 Data/all_samples_Data_Repo_MiSeq2_10c_complete.tsv | while IFS= read -r sample; do
	if grep -q "$sample" Results/Fq_metadata/rawfilelist.tsv 
	then
		echo "${sample} present in raw files"
	else
		echo "$sample" >> Results/missing_file.tsv
	fi; 
done
