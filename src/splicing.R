# Install and load the GenomicRanges package
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges", force = TRUE)
library(GenomicRanges)

# Function to identify exons and introns
find_exons_introns <- function(nucleotide_seq) {
    # Count the number of nucleotides in the sequence
    seq_length <- nchar(nucleotide_seq)

    # Create a vector to store exon and intron positions
    positions <- character(seq_length)
    positions <- rep("I", seq_length) # Initialize all positions as introns

    # Identify exons
    for (i in 1:seq_length) {
    # If the current nucleotide is not an intron, it's an exon
    if (substr(nucleotide_seq, i, i) != "I") {
        positions[i] <- "E"
    }
    }

    # Create GRanges objects for exons and introns
    exons_gr <- GRanges(seqnames = "seq1", ranges = IRanges(start = 1, end = seq_length))[positions == "E"]
    introns_gr <- GRanges(seqnames = "seq1", ranges = IRanges(start = 1, end = seq_length))[positions == "I"]

    # Print exons and introns
    print("Exons:")
    print(exons_gr)
    print("Introns:")
    print(introns_gr)
}

# Example usage of the function
nucleotide_sequence <- "AUGGUAUAGCUAAUGGU"
find_exons_introns(nucleotide_sequence)
