#!/bin/R env

########################################
# LOAD LIBRARIES AND PREPARE VARIABLES #
########################################

if (!require("argparse")) {
  install.packages("argparse", repos="http://cran.rstudio.com/")
  library("argparse")
}

if (!require("stringdist")) {
  install.packages("stringdist", repos="http://cran.rstudio.com/")
  library("stringdist")
}

# Custom filtering, denoising parameters (if not default) can be provided as a separate config file?
parser <- ArgumentParser()
parser$add_argument("-p", "--path_to_fastq", help="Path to merged fastq file (required)")
parser$add_argument("-d", "--dir", help="Working directory path for writing all bbmerge contamination output files")
parser$add_argument("-b", "--barcodes", help="The path to a csv file with the sample_id,Forward,Reverse, where Forward and Reverse are columns with the barcodes for the sample")
args <- parser$parse_args()

# Universal parameters
fastq_file <- args$path_to_fastq
work_dir <- args$dir
path_to_flist <- args$barcodes

print(paste0("Processing file: ", fastq_file))
write(paste0("Processing file: ", fastq_file), stderr())
########################################
#             BARCODE REPORT           #
########################################

source(file.path("/Code", "matching_functions.R"))
#source(paste0(file.path(dirname(dirname(work_dir)), "Code", "matching_functions.R")))
barcodes = read.csv(path_to_flist, sep = ",", header = TRUE)
dist = 2

# Open the Fastq file
con <- file(fastq_file, "r")
sample <- sub(".*/([^/]+)_merged\\.fastq", "\\1", fastq_file)
  
seq_names = c()
forward_barcodes = c()
reverse_barcodes = c()
forward_distances = c()
reverse_distances = c()
five_prime_ends = c()
three_prime_ends = c()
match_categories = c()
insert_sizes = c()
print(sample)
  
while (length(line <- readLines(con, n = 4)) > 0) {
  seq_name <- line[1]
  sequence <- line[2]

  match_tmp = match_fun(sample, sequence, barcodes, dist)
    
  seq_names = c(seq_names, seq_name)
  forward_barcodes = c(forward_barcodes, match_tmp[[1]])
  reverse_barcodes = c(reverse_barcodes, match_tmp[[2]])
  forward_distances = c(forward_distances, match_tmp[[3]])
  reverse_distances = c(reverse_distances, match_tmp[[4]])
  five_prime_ends = c(five_prime_ends, match_tmp[[5]])
  three_prime_ends = c(three_prime_ends, match_tmp[[6]])
  match_categories = c(match_categories, match_tmp[[7]])
  insert_sizes = c(insert_sizes, match_tmp[[8]])
    
}
  
close(con)
  
df = data.frame(sequence = seq_names,
                 forward_barcodes = forward_barcodes,
                 reverse_barcodes = reverse_barcodes, 
                 forward_distances = forward_distances,
                 reverse_distances = reverse_distances,
                 five_prime_ends = five_prime_ends,
                 three_prime_ends = three_prime_ends,
                 match_status = match_categories,
                 insert_sizes = insert_sizes)
              
output_filename_o = file.path(dirname(dirname(work_dir)), "Report", "Merge", paste0(sample, "_final.tsv"))
write.table(df, file=output_filename_o, quote = FALSE, sep = "\t", row.names = FALSE)
print(paste0("BBmerge contamination report written to ", output_filename_o))
