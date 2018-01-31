# DNAreadsToBytes
Decodes digital (byte) information encoded in short DNA sequences (fasta) according to Goldman et al. (2013) and saves it as files.

DNAreadsToByte is a simple R script that takes a fasta file of short 117-pb long reads that are sequenced from artificially produced DNA fragments used to encode computer files, according to Goldman et al. Nature. (2013) doi:  10.1038/nature11875. Decoding follows the specification 2.0 as detailed in https://www.ebi.ac.uk/sites/ebi.ac.uk/files/groups/goldman/file2features_2.0.pdf

The 117-bp requirement represents the length of the artifical fragments, so the user should ensure that all sequences of the fasta are exactly this long, or else the script may not behave properly.

The script needs the path of the fasta file, and (optionally) the directory where the user wants the decoded files to be written. File names use the prefix "file" followed by a number from 0 to 8, with no extension. The user will have to figure out the types of the generated files. 

As the complete fasta is imported in the R environment, the user should ensure that their computer has enough free RAM for the script to run, that is, at least 2-3 times the size of the fasta file.

The script has been tested in R 3.3.2 and requires the following packages:
data.table
stringi
Biostrings
matrixStats


