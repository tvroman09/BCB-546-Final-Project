library(Biostrings)

fasta <- readDNAStringSet("NPPF2.fasta") #read in FASTA file


library(msa)

align <- msa(fasta) #multiple sequence alignment of all seqs from FASTA file
align <- DNAStringSet(align) #store the alignments as DNAString objects

writeXStringSet(align, "alignment.fasta") #writing the alignment to a new file

library(ape)
library(phangorn)

# Read the alignment file 
alignment <- read.phyDat("alignment.fasta", format = "fasta") 

# Build a distance matrix from the alignment
dist_matrix <- dist.hamming(alignment) 

# Build the phylogenetic tree using neighbor-joining method
phy_tree <- nj(dist_matrix) 

# Visualize the phylogenetic tree
plot(phy_tree) 
