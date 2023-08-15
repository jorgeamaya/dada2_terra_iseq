#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
data_dir = args[1]
out_dir =args[2]

data_dir = "/Users/jorgeamaya/Desktop/Broad_Test/amplicon_decontamination_pipeline/Report/Merge/" 
out_dir = "/Users/jorgeamaya/Desktop/Broad_Test/amplicon_decontamination_pipeline/Report/"

if (!require("viridis")) {
  install.packages("viridis", repos="http://cran.rstudio.com/")
  library("viridis")
}

###################################
###BB MERGE PERFORMANCE REPORT ####
###################################

mergedata <- read.csv(file.path(data_dir, "bbmergefields.tsv"), sep = "\t", header = TRUE)
subset = read.csv(file.path(data_dir, "subset.tsv"), sep = "\t", header = TRUE)[[1]]

mergedata = mergedata[mergedata$SampleID %in% subset,]

pdf(file.path(out_dir, "BBmerge_performance_subset_report.pdf"), width = 12)
barplot(mergedata[,2] ~ mergedata[,1], las = 2, cex.names = 1)
dev.off()
