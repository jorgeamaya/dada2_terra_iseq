#!/bin/r env

########################################
# LOAD LIBRARIES AND PREPARE VARIABLES #
########################################

if (!require("seqinr")) {
  install.packages("seqinr", repos="http://cran.rstudio.com/")
  library("seqinr")
}

print('LOADED THE FIRST LIBRARY')
#if (!require("data.table")) {
#  install.packages("data.table", repos="http://cran.rstudio.com/")
#  library("data.table")
#}
#if (!require("argparse")) {
#  install.packages("argparse", repos="http://cran.rstudio.com/")
#  library("argparse")
#}
#if (!require("Biostrings")) {
#  install.packages("Biostrings", repos="http://cran.rstudio.com/")
#  library("Biostrings")
#}
#if (!require("parallel")) {
#  install.packages("parallel", repos="http://cran.rstudio.com/")
#  library("parallel")
#}
#if (!require("doMC")) {
#  install.packages("doMC", repos="http://cran.rstudio.com/")
#  library("doMC")
#}
#
#absolute = function(correctedASV, tar) {
#  return(nchar(correctedASV) != nchar(tar))
#}
#per_size = function(correctedASV, tar) {
#  return((abs(nchar(correctedASV)-nchar(tar))/nchar(tar)) > 0.1)
#}
#levenshtein = function(correctedASV, tar) {
#  positions <- list(first = gregexpr("N", correctedASV)[[1]][1],
#                    last = tail(gregexpr("N", correctedASV)[[1]], 1))
#  correctedASV_three_prime = substr(correctedASV, 1, positions$first-1)
#  correctedASV_five_prime = substr(correctedASV, positions$last+1, nchar(correctedASV))
#  ref_three_prime = substr(tar, 1, positions$first-1)
#  ref_five_prime = substr(tar, positions$last+1, nchar(tar))
#  
#  return((stringDist(c(paste0(correctedASV_three_prime, correctedASV_five_prime),
#                      paste0(ref_three_prime, ref_five_prime)),
#                    method="levenshtein")) > 9)
#}
#
#parser <- ArgumentParser()
#parser$add_argument("-s", "--seqtab", 
#                    help="Path to input")
#parser$add_argument("-ref", "--reference",
#                    help="Path to reference fasta sequences")
#parser$add_argument("-dist", "--distance",
#                    help=paste("Distance method to accept the corrected ASV. It can be absolute, percentage, levenshtein, or ignore",
#                               "absolute: Corrected ASV discarded if its size is different to the reference.",
#                               "percentage: Corrected ASV discarded if its size is 10 larger or smaller than the reference.",
#                               "levenshtein: Corrected ASV discarded if its size is 9 bp larger or smaller than the reference.",
#                               "ignore: Corrected ASV never discarded.", sep = " "))
#parser$add_argument("-o", "--output",
#                    help="Path to output for corrected ASV list")
#
#args <- parser$parse_args()
#path_to_refseq <- args$reference
#seqfile <- args$seqtab
#output <- args$output
#dist = args$dist
#
#if (dist == 'ignore') {
#  print("Caution: All Corrected ASV will be included in the output. Some ASVs in your table may be incorrect constructs.")
#}
#
##seqtab = '/Desktop/Broad_Test/amplicon_decontamination_pipeline/Results/DADA2_NOP/seqtab.tsv'
##reference= '~/Desktop/Broad_Test/amplicon_decontamination_pipeline/Data/pf3d7_ref_updated_v4.fasta'
#
#if (file.exists(path_to_refseq)) {
#  ref <- toupper(sapply(read.fasta(path_to_refseq),c2s))
#} else {
#  stop("Reference file not found!")
#}
#
#if (file.exists(seqfile)) {
#  seqtab <- as.matrix(fread(seqfile), rownames=1)
#} else {
#  stop(paste("ASV sequence table file",seqtab,"not found!"))
#}
#
#########################################
##           PROCESS ASVs               #
#########################################
#
##Generate the substitution matrix
#sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
#seqs <- as.character(colnames(seqtab))
#
#registerDoMC(detectCores())
#df <- foreach(i=1:length(seqs), .combine = "rbind") %dopar% {
#  map <- pairwiseAlignment(ref, seqs[i], substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = TRUE)
#  tar = ref[which.max(map)]
#  seq <- strsplit(seqs[i],"NNNNNNNNNN")[[1]]
#  aln <- pairwiseAlignment(seq[1:2], tar, substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE, type = 'overlap')
#  con <- compareStrings(consensusString(aln[1]), consensusString(aln[2]))
#  overlap <- unlist(gregexpr("[[:alpha:]]", con))
#  
#  if (overlap == -1) {
#    N = (nchar(seq[1])+nchar(seq[2])) - nchar(tar)
#    stkN <- paste0(rep('N', abs(N)), collapse = '')
#    correctedASV <- paste0(seq[1], stkN, seq[2])
#  } else {
#    N = length(overlap)
#    correctedASV <- paste0(seq[1], substr(seq[2], (N+1), nchar(seq[2])))
#  }
# 
#  if(dist == 'absolute') {
#    if (absolute(correctedASV, tar)) {
#      N = NA
#      correctedASV = NA
#    }
#  } else if (dist == 'levenshtein') {
#    if (levenshtein(correctedASV, tar)) {
#      N = NA
#      correctedASV = NA
#    }
#  } else if (dist == 'percentage') {
#    if (absolute(correctedASV, tar)) {
#      N = NA
#      correctedASV = NA
#    }
#  }
#  
#  data.frame(target = names(tar),
#             ASV = seqs[i],
#             correctedASV = correctedASV,
#             overlap = N)
#}
#
#print(df)
#write.table(df, file = output, sep = "\t", quote = FALSE, row.names = FALSE)
#seqfile_corrected <- paste0(dirname(seqfile), "/seqtab_corrected.tsv")
#colnames(seqtab) <- as.character(df$correctedASV)
#	#seqtab = seqtab[,which(colnames(seqtab) == NA)]
#print(seqtab)
#	#rownames(seqtab) <- as.character(rownames(df))
#write.table(seqtab, file = seqfile_corrected, sep = "\t", quote = FALSE, row.names = TRUE)
