#!/bin/r env

########################################
# LOAD LIBRARIES AND PREPARE VARIABLES #
########################################

if (!require("argparse")) {
  install.packages("argparse", repos="http://cran.rstudio.com/")
  library("argparse")
}
if (!require("limma")) {
  install.packages("limma", repos="http://cran.rstudio.com/")
  library("limma")
}
if (!require("data.table")) {
  install.packages("data.table", repos="http://cran.rstudio.com/")
  library("data.table")
}
if (!require("stringr")) {
  install.packages("stringr", repos="http://cran.rstudio.com/")
  library("stringr")
}
if (!require("seqinr")) {
  install.packages("seqinr", repos="http://cran.rstudio.com/")
  library("seqinr")
}
if (!require("parallel")) {
  install.packages("parallel", repos="http://cran.rstudio.com/")
  library("parallel")
}
if (!require("doMC")) {
  install.packages("doMC", repos="http://cran.rstudio.com/")
  library("doMC")
}
if (!require("Biostrings")) {
  install.packages("Biostrings", repos="http://cran.rstudio.com/")
  library("Biostrings")
}

# Pairwise Alignment
seq_align <- function(seqs_df, path_to_ref, overlap = TRUE, parallel = TRUE) {
  align_df <- data.frame()
  # LOAD REFERENCE
  if (file.exists(path_to_ref)) {
    ref <- read.fasta(path_to_ref)
    ref_str <- toupper(sapply(ref, c2s))
  } else {
    stop(paste("File", path_to_ref, "not found!"))
  }

  sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
  if (!overlap) {
    ###### TO DO ######
    #Deal with non-ovlerpapping iseq reads. For the moment, this function seems to work with uncorrected ASVs, i.e.,
    #ASVs that have 10 Ns by default. This function will need working. It's throwing the following error.
    #strsplit2(x = 'AAAXXXXXXXXXXAA', split = split)
    #Error in as.vector(x, "character") : 
    #cannot coerce type 'closure' to vector of type 'character'
    ###################
    
    split <- paste0(rep("N", 10), collapse = "")
    seq_split <- strsplit2(x = seqs_df[,1], split = split)
    seq_all <- data.frame(sequence=paste0(seq_split[,1],seq_split[,2]), hapid = seqs_df[,2])
  } else {
    seq_all <- data.frame(sequence=seqs_df[,1], hapid = seqs_df[,2])
  }

  if (parallel) {
    registerDoMC(detectCores())
    align_df <- foreach(seq_1=1:length(seq_all$sequence), .combine = "rbind") %dopar% {
      #Find out the best refrence sequence to which each ASV aligns best. 
      aln <- pairwiseAlignment(ref_str, seq_all$sequence[seq_1], substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
      num <- which.max(score(aln))
      patt <- c(alignedPattern(aln[num]), alignedSubject(aln[num])) #Pattern is the reference and Subject is the ASV.
      #Obtain the generalized Levenshtein distance, i.e., the minimum number of inversion, 
      #deletetions, and substitution to transfrom one string into another
      dist <- adist(as.character(patt)[1],as.character(patt)[2])
      ind <- sum(str_count(as.character(patt),"-"))
      data.frame(hapid = seq_all$hapid[seq_1], #Haplotype id as ASV1, ASV2, ASV3,
                 hapseq = as.character(patt)[2], #Haplotype/ASV 
                 refseq = as.character(patt)[1], #Reference Sequence
                 refid = names(patt)[1], #Reference Sequence ID
                 aln_score = score(aln[num]), #Score of the alignment (the best alignment) 
                 snv_dist = as.integer(c(dist - ind)), #How much of the Levenshtein distances are subtitutions.
                 indel_dist = as.integer(ind)) #Number of indels
    }
  } else {
    for (seq_1 in 1:length(seq_all$sequence)) {
      aln <- pairwiseAlignment(ref_str, seq_all$sequence[seq_1], substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
      num <- which.max(score(aln))
      patt <- c(alignedPattern(aln[num]), alignedSubject(aln[num]))
      dist <- adist(as.character(patt)[1],as.character(patt)[2])
      ind <- sum(str_count(as.character(patt),"-"))
      df <- data.frame(hapid = seq_all$hapid[seq_1], 
                       hapseq = as.character(patt)[2], 
                       refseq = as.character(patt)[1], 
                       refid = names(patt)[1], 
                       aln_score = score(aln[num]), 
                       snv_dist = as.integer(c(dist - ind)),
                       indel_dist = as.integer(ind))
      align_df <- rbind(align_df,df)
    }
  }
  return(align_df)
}

parser <- ArgumentParser()
#Minimum arguments
parser$add_argument("-s", "--seqtab", help="Path to input")
parser$add_argument("-b","--bimera", help="ASV File with identifed bimeras")
parser$add_argument("-snv", "--snv_filter", help="Path to file for filtering ASVs based on edit distance")
parser$add_argument("--indel_filter", help="Specify proportion of ASV length (between 0 and 1) to target length for filtering based on indels")
parser$add_argument("-o", "--output", help="Path to output file")
parser$add_argument("--fasta", action='store_true', help="Write ASV sequences separately into fasta file")
parser$add_argument("--parallel", action='store_true', help="Enable parallel processing")
#Variable arguments
parser$add_argument("-noref", "--no_reference", action='store_true', help="specify no reference provided")
parser$add_argument("-ref", "--reference", help="Path to reference fasta sequences")
parser$add_argument("--strain", default="3D7", help="Name of Specific strain to map to. Defaults to 3D7")
parser$add_argument("-ref2", "--reference2", help="Path to reference2 fasta sequences")
parser$add_argument("--strain2", help="Name of second strain if mapping to 2 different strains")

args <- parser$parse_args()

########################################
#           PROCESS ASVs               #
########################################

# PostProc Begin
output <- args$output
seqfile <- args$seqtab
path_to_refseq <- args$reference
strains <- args$strain
parallel = args$parallel
path_to_refseq2 = args$reference2
strain2 = args$strain2
no_reference = args$no_reference
filter_file = args$snv_filter
indel_filter = args$indel_filter
bimera = args$bimera
fasta = args$fasta
output = args$output

#output = '/Users/jorgeamaya/Desktop/dada2_terra_iseq/Results/PostProc_DADA2/ASVTable.txt'
#seqfile = '/Users/jorgeamaya/Desktop/dada2_terra_iseq/Results/seqtab_mixed.tsv'
#path_to_refseq = 'Data/pf3d7_ref_updated_v4.fasta'
#strains = '3D7'
#strain2 = 'DD2'
#parallel = FALSE
#path_to_refseq2 = 'Data/pfdd2_ref_updated_v3.fasta'
#no_reference = FALSE
#filter_file = 'Data/snv_filters.txt'
#indel_filter = 0.895
#bimera = '/Users/jorgeamaya/Desktop/dada2_terra_iseq/Results/ASVBimeras.txt'

#LOAD THE REFERENCES AND THE STRAINS
if (!no_reference) {
  print('Running postProc with reference.')
  if (is.null(args$reference) && is.null(args$strain)) {
    stop("Reference genome (--reference) name of target strain (--strain) missing.")
  }

  if (!is.null(args$reference2) || !is.null(args$strain2)) {
    print('Running postProc with second reference.')
    if (is.null(args$reference2) && is.null(args$strain2)) {
      stop("Reference genome (--reference2) name of target strain (--strain) missing.")
    } else {
      path_to_refseq <- c(path_to_refseq, path_to_refseq2)
      strains <- c(strains, strain2)
    }
  }
}
print("#Produce the asvdf and seqs_df")
#Produce the asvdf and seqs_df
if (file.exists(seqfile)) {
	seqtab <- as.matrix(fread(seqfile), rownames=1)
	# Identify columns with names shorter than 4 letters
	cols_to_remove <- colnames(seqtab)[nchar(colnames(seqtab)) < 10]
	# Remove identified columns
	seqtab <- seqtab[, !(colnames(seqtab) %in% cols_to_remove)]

  seqs <- colnames(seqtab)
	nsample=nrow(seqtab)
	hapid <- paste0("ASV",1:length(seqs))
	# DataFrame for aligning to truth set
	seqs_df <- data.frame(sequence = seqs, hapid = hapid)
	# Change colnames of ASV from sequences to ASV ids
	seqtab_haps <- seqtab
	colnames(seqtab_haps) <- hapid
	#seqtab and seqtab_haps are the same matrix except that the colnames have been changed
	#from the sequence of the ASVs to happlotype ids of the form ASV1, ASV2, ASV3, ...
	
	## ASV summary table
	total_reads <- apply(seqtab_haps,2,sum) #Sum column wise to obtain the total number of occurrences of each ASV across samples
	total_samples <- apply(seqtab_haps,2,function(x) sum(x != 0)) #Obtain the total number of samples/rows in which the ASV is present.
	asvdf <- data.frame(hapid = hapid, #hapid in the ASV1, ASV2, ASV3, ... form.
	                    haplength = nchar(seqs), #Length of the ASV.
	                    total_reads = total_reads, #Total number of occurrences of each ASV across samples
	                    total_samples = total_samples, #Total number of samples/rows in which the ASV is present.
	                    strain = "N") #Strain column. Dummy column for the moment. See Map True Set commands 
	asvdf$hapid <- as.character(asvdf$hapid) #This needs to be done to prevent an invalid factor level, NA generated error
	asvdf$strain <- as.character(asvdf$strain) 
} else {
  stop(paste("ASV sequence table file", seqtab, "not found!"))
}

print("no reference")
if (!no_reference) {
	for (p in 1:length(path_to_refseq)) {
	  # Alignment with RefSet
	  align_df <- seq_align(seqs_df = seqs_df, path_to_ref = path_to_refseq[p], overlap = TRUE, parallel = parallel) #This overlap will have to be modifiedto accomodate mixed_reads
	  #align_df$refid <- as.character(align_df$refid) #Eliminate
	  #align_df$hapid <- as.character(align_df$hapid) #Eliminate
	  
	  ## Map True Set onto ASV summary table based on exact and inexact matches to true set
	  df <- align_df[,c(1,4,6,7)]
	  colnames(df) <- c("hapid", paste0("refid_", strains[p]), paste0("snv_dist_from_", strains[p]), paste0("indel_dist_from_", strains[p]))
	  asvdf <- merge(asvdf, df, by = "hapid", sort = FALSE)
	  #Strain = The name if the strain in which the ASV is found. Snv dist and indel dist must equal 0.
	  snv_dist_l = as.logical(as.numeric(align_df$snv_dist) == 0)
	  indel_dist_l = as.logical(as.numeric(align_df$indel_dist) == 0)
	  asvdf$strain[snv_dist_l & indel_dist_l] <- as.character(strains[p]) #If a strain[p] has alredy being declared, this step will overwrite it.
	}

  #Filter against SNV filter
  if (!is.null(args$snv_filter)) {
    if (file.exists(filter_file)) {
      VariantCounts <- fread(filter_file)
  		asvdf$snv_filter <- NA
    	# For Sliding edit distance and length based filter
  		for (i in 1:nrow(asvdf)) {
  		  refid <- asvdf[i,paste0("refid_", strains[1])]
  		  #hapdist = the snv distance of the haplotype
      	hapdist <- asvdf[i,paste0("snv_dist_from_", strains[1])]
      	#tardist = the variant counts Ref Sequence that matches best the ASV
  			tardist <- VariantCounts[c(VariantCounts$id == refid),2]
   			if (hapdist <= tardist) {
    			sfil <- "PASS"
  			} else {
    			sfil <- "FAIL"
   			}
  			asvdf$snv_filter[i] <- sfil #The haplotype is way too distant from the target ASV, even after allowing the SNV filter changes.
    	}
  	} else {
  		warning(paste("File",filter_file,"not found!. Skipping SNV based filtering.."))
  	}
	} else {
  	print("SNV filter file not given. Skipping SNV based filtering..")
	}
  print("After SNV filter")
  
  #Filter against indel criteria
  if (!is.null(args$indel_filter)) {
  	if (file.exists(indel_filter)) {
  		Indelcutoff <- fread(indel_filter)
		} else {
   		Indelcutoff <- as.numeric(indel_filter)
  	}
  	if (is.na(Indelcutoff)) {
  		warning("INDEL based filter threshold argument not valid. Using default proportion of +-90%..")
  		Indelcutoff <- 0.90
  	}
    
		asvdf$indel_filter <- NA
		refseq <- toupper(sapply(read.fasta(path_to_refseq[1]),c2s))
  	for (i in 1:nrow(asvdf)) {
  		haplen <- nchar(as.character(seqs_df$sequence[i]))
  		tarlen <- nchar(refseq[refid])
  		lprop <- as.numeric(haplen/tarlen)
  		
  		refid <- asvdf[i,paste0("refid_",strains[1])]
  		########This if statement processes the idnel cutoff when provided with table.
  		if (length(Indelcutoff) != 1) {
  			ind <- Indelcutoff[c(Indelcutoff[,1] == refid),2]
  		} else {
  			ind <- Indelcutoff
  		}
  		
  		if (lprop < ind | lprop > as.numeric(1/ind)) { #(the asv is ind*100 perc shorter than the reference | the asv is ind*100 perc larger than the reference. Notice that this last case is less precise because the 1/ind may be periodic  )
  			ifil <- "FAIL"
  		} else {
  			ifil <- "PASS"
  		}
  		asvdf$indel_filter[i] <- ifil #The sample ASV is way shorter or way longer than the reference ASV
		}
	} else {
  		print("INDEL based filter threshold not given. Skipping length filtering..")
	}
  print("After Indel Filter")
}
print("bimera")  
#Add the bimera information by matching the bimera table to asvdf
if (!is.null(args$bimera)) {
  if (file.exists(bimera)) {
    bimeras <- fread(bimera)
    bimeras_tmp = bimeras[,c(4,5)]
    bimeras = bimeras[,c(2,3)]
    bimeras = rbind(bimeras, bimeras_tmp) 
    seqs_df <- merge(seqs_df,bimeras, by = "sequence", all = TRUE, sort = FALSE)
    asvdf <- merge(asvdf,seqs_df[,2:3], by = "hapid", all = TRUE, sort = FALSE)
  } else {
    warning(paste("File", bimera, "not found. Skipping bimera flag.."))
  }
} else {
  print("Bimeric ASV file not given. Skipping bimera flag..")
}

asvdf = asvdf[!is.na(asvdf$hapid),]
#A table with the columns
#"hapid", "haplength", "total_reads", "total_samples", "strain", "refid_3D7", "snv_dist_from_3D7", 
#"indel_dist_from_3D7", "refid_DD2", "snv_dist_from_DD2", "indel_dist_from_DD2", "snv_filter", "indel_filter" 
write.table(asvdf, file = output, sep = "\t", quote = FALSE, row.names = FALSE)

#A fasta fiel with the ASV with identifiers as ASV1, ASV2, ASV3, ...
if (fasta) {
  write.fasta(sapply(seqs, s2c), names = hapid, file.out = paste0(dirname(output),"/ASVSeqs.fasta"), nbchar = 600)
}
