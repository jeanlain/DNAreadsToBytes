# DNAreadsToBytes
Decodes digital (byte) information encoded in short DNA sequences (fasta) according to Goldman et al. (2013) and saves it as files.

DNAreadsToByte.R is a simple R script that takes a fasta file of  117-bp-long reads that are sequenced from artificially produced DNA fragments used to encode computer files, according to Goldman et al. Nature. (2013) doi:  10.1038/nature11875. Decoding follows the specification 2.0 as detailed in https://www.ebi.ac.uk/sites/ebi.ac.uk/files/groups/goldman/file2features_2.0.pdf

The 117-bp requirement represents the length of the artifical fragments, so the user should ensure that all sequences of the fasta are exactly this long, or else the script may not behave properly.

usage:

Rscript   DNAreadsToBytes.R   path_to_reads.fasta   folder_where_you_want_decoded_files   path_to_huffman_code_table_file

If the last argument is not provided, the script will look for the file "View_huff3.cd.new.correct" in the directory specified as the second argument.

Decoded file names will use the prefix "file" followed by a number from 0 to 8, with no extension. 

As the complete fasta is imported in the R environment, users should ensure that their computer has enough free RAM for the script to run, that is, at least 3 times the size of the fasta file.

DNAreadsToBytes.R has been tested in R 3.3.2 and requires the following packages to be installed:

data.table
stringi
Biostrings
matrixStats


