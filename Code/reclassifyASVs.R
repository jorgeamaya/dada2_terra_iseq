#!/bin/r env

########################################
# LOAD LIBRARIES AND PREPARE VARIABLES #
########################################

if (!require("argparse")) {
  install.packages("argparse", repos="http://cran.rstudio.com/")
  library("argparse")
}
if (!require("stringr")) {
  install.packages("stringr", repos="http://cran.rstudio.com/")
  library("stringr")
}

count_Ns <- function(seqtab) {
  num_Ns <- sapply(seqtab, function(x) {
    Ns = unlist(gregexpr("N", x, fixed = TRUE))
    if (length(Ns) == 1) {
      if (Ns == -1) {
        return(0)
      } else {
        return(1)
      }
    } else {
      return(sum(lengths(Ns)))
    }
  })
  return(num_Ns)
}

#parser <- ArgumentParser()
#parser$add_argument("-s", "--seqtab", 
#                    help="Path to input")
#parser$add_argument("-o", "--output",
#                    help="Path to output for corrected ASV list")

#seqfile <- args$seqtab
#output <- args$output

#########################################
###           LOAD DATA               ###
#########################################

seqfile = '/Users/jorgeamaya/Desktop/Broad_Test/amplicon_decontamination_pipeline/Results/seqtab_iseq.tsv'

if (file.exists(seqfile)) {
  seqtab <- as.matrix(fread(seqfile), rownames=1)
} else {
  stop(paste("ASV sequence table file",seqtab,"not found!"))
}

#Correct the first NA. The first NA has no .[d+] append at the end. The subsequent algorithms will not detect it.
colnames(seqtab)[grepl("^V\\d+$", colnames(seqtab))] = 'NA.0'

df = data.frame(id=1:ncol(seqtab),
                ASVs = colnames(seqtab),
                NAs = grepl("^NA\\.\\d+$", colnames(seqtab)),
                size = nchar(colnames(seqtab)),
                N_number = count_Ns(colnames(seqtab))
                ) #Notice that simply detecting NA may match some sequences because N is used to pad the ASVs.

df$size[df$NAs] = rep(NA, sum(df$NAs)) 
df$N_number[df$NAs] = rep(NA, sum(df$NAs)) 
df$groups = rep(NA, nrow(df))

dfNs = df[df$NAs,]
df = df[!df$NAs,]

#Order sequences by size
df = df[order(df[,4], df[,5]), ]

levenshtein_matrix = adist(df$ASVs)

num_diff_Ns <- matrix(0, nrow = length(df$ASVs), ncol = length(df$ASVs))
for (i in 1:length(df$ASVs)) {
  for (j in 1:length(df$ASVs)) {
    if (i != j) {
      num_diff_Ns[i, j] <- abs(sum(unlist(gregexpr("N", df$ASVs[i])) > 0) - sum(unlist(gregexpr("N", df$ASVs[j])) > 0))
    }
  }
}

#Find the amount of Levensthein difference not explained by indels
#In other words, if the sequences have MORE DIFFERENCES than their absolute number of Ns differences,
#then they are different ASVs and must be preserved as independent. On the contrary, if all the difference
#they have is explained by different number of N, then they are the same ASV.
Ns_matrix = levenshtein_matrix - num_diff_Ns
TrueFalse_Ns_matrix = (Ns_matrix == 0)
colnames(TrueFalse_Ns_matrix) = df$ASVs
rownames(TrueFalse_Ns_matrix) = df$ASVs

#Make the lower diagonal FALSE to avoid duplication when melting the matrix 
TrueFalse_Ns_matrix[lower.tri(TrueFalse_Ns_matrix)] = FALSE
#Make the diagonal FALSE to remove self comparisons
diag(TrueFalse_Ns_matrix) = FALSE

#Melt the Matrix
TrueFalse_Ns_long = melt(TrueFalse_Ns_matrix)

#Find the pairs that have different numbers of N but are the same sequence
pairs = TrueFalse_Ns_long[TrueFalse_Ns_long$value,] 

#Create a NULL vector for the groups and initialize the first pair
pairs$group = NULL
pairs$group[1] = 1

#Some pairs may be in the table more than once. Consider,
#1,ASV1,ASV2
#2,ASV2,ASV3
#ASV1 == ASV2 and ASV2 and ASV3; therefore, ASV1 == ASV3.
#They all belong to the same group.

#This for loop compares all possible pairs and to find this groupings.
for (tar_pair in 2:nrow(pairs)) {
  found = FALSE
  for (ref_pair in 1:c(tar_pair -1)) {
    if(pairs[ref_pair,1] == pairs[tar_pair,1] |
       pairs[ref_pair,1] == pairs[tar_pair,2] | 
       pairs[ref_pair,2] == pairs[tar_pair,1] |
       pairs[ref_pair,2] == pairs[tar_pair,2] ) {
      pairs$group[tar_pair] = pairs$group[ref_pair]   
      found = TRUE
      break
    }
  }
  if(!found) {
    print(unique(pairs$group))
    pairs$group[tar_pair] = max(unique(pairs$group))+1
  }
}

cppairs = data.frame(ASV=c(pairs$X1, pairs$X2),
                     group=c(pairs$group, pairs$group))

cppairs = cppairs[!duplicated(cppairs), ]

cppairs$Ns = sapply(cppairs$ASV, function(x) length(unlist(gregexpr("N", x, fixed = TRUE))))

seqtab <- as.matrix(fread(seqfile), rownames=1)
for (group in 1:35) { #unique(cppairs$group)) {
  #Find the ASV with the least number of Ns
  group = 36
  ASVs = cppairs[cppairs$group == group,]
  ASVleastN = ASVs$ASV[which(ASVs$Ns == min(ASVs$Ns))]
  colnames(seqtab)[which(colnames(seqtab) %in% ASVs$ASV)] = as.character(rep(ASVleastN, sum(colnames(seqtab) %in% ASVs$ASV)))
}

collapsed_seqtab <- as.data.frame(t(rowsum(t(seqtab), group = colnames(seqtab))))

seqfile_corrected <- paste0(dirname(seqfile), "/seqtab_collapsed.tsv")
write.table(collapsed_seqtab, file = seqfile_corrected, sep = "\t", quote = FALSE, row.names = TRUE)
