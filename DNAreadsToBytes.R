

require(stringi)
require(data.table)
require(Biostrings)
require(matrixStats)
options("stringsAsFactors" = F)

paths <- commandArgs(T)

# path to fasta file
input <- paths[1]

# output directory
outputFolder <- paths[2]

if (length(paths) < 3) {
  huffmanCode <- gsub("//", "/", stri_c(outputFolder, "/View_huff3.cd.new.correct"), fixed = T)
} else {
  huffmanCode <- paths[3]
}

paths <- c(input, outputFolder, huffmanCode)
missing <- !file.exists(paths)

if (any(missing)) stop(paste("missing file(s):", paste(paths[missing], collapse = ", ")))



# defining "basic" functions" used in the script-----------------------------------------

baseTypes <- c("A", "C", "G", "T")


# matrix used to convert successive bases to trits, used below
conv <- matrix(c(NA, 0:2, 2L, NA, 0:1, 1:2, NA, 0L, 0:2, NA), 4, byrow = T)
colnames(conv) <- baseTypes
rownames(conv) <- baseTypes


# converts a matrix of DNA bases (sequences as rows) to integer trits according to Goldman et al's scheme
toTrits <- function(mat) {

  # in case mat is a vector, we convert it to a 1-row matrix (otherwise we would have dimension issues)
  if (!is.matrix(mat)) mat <- rbind(mat)
  vapply(2:ncol(mat), function(col) conv[cbind(mat[, col - 1], mat[, col])], integer(nrow(mat)))
}


# splits a string vector into a matrix with one character per cell. Words are
# rows. The number of columns is the length of the longest word. We use as
# substring function (stri_sub()), which appears to be faster than a split
stringToMatrix <- function(vector) {
  nchars <- stri_length(vector)
  m <- max(nchars)
  mat <- matrix(NA, nrow = length(vector), ncol = m, byrow = T)
  for (i in 1:m) {
    mat[, i] <- stri_sub(vector, i, i)
  }
  mat
}


# reverse complements a character vector representing DNA sequences. Ambiguities are not treated
revCom <- function(seqs) {
  stri_reverse(chartr("acgtACGT", "tgcaTGCA", seqs))
}



# converts numbers of a given base to base 10, assuming digits are in separate
# columns of a matrix, so that rows represent numbers read from left to right
toBase10 <- function(mat, base = 3) {
  if (!is.matrix(mat)) mat <- rbind(mat)
  res <- 0
  ncol <- ncol(mat)
  for (col in ncol:1) {
    res <- res + mat[, col] * base^(ncol - col)
  }
  res
}


# tells the occurrence of the element found, for all positions of a vector (e.g.
# 1 if it's the first time it appears, 2 if the 2nd time, etc)
occurrences <- function(vect) {
  dt <- data.table(vect, pos = 1:length(vect))
  counts <- dt[, .(occ = 1:.N, pos), by = vect]
  res <- integer(length(vect))
  res[counts$pos] <- counts$occ
  res
}


# returns the first position (column) of each element of "values" in each
# corresponding row of matrix "mat" (recycling may work, as values can have 1
# element)
rowMatches <- function(values, mat) {
  res <- which(mat == values, arr.ind = T)
  res[match(1:nrow(mat), res[, 1]), 2]
}


# returns a matrix in which colums represent the cumulated sums of the colums of mat, starting from the left
cumsumMatrix <- function(mat) {
  mat <- cbind(mat, -rowSums(mat))
  cs <- matrix(cumsum(t(mat)), nrow(mat), byrow = T)
  cs[, -ncol(cs)]
}


# splits a character vector into columns (like excel). A column range can be
# specified as an integer vector and "cells" are filled with NAs if needed. By
# default, it generate as many columns as necessary given the character(s) used
# to split the vector elements.
splitToColumns <- function(vector, split, columns) {
  vector <- as.character(vector)
  res <- stri_split_fixed(vector, split, omit_empty = NA, simplify = T)
  if (missing(columns)) {
    columns <- 1:ncol(res)
  }
  if (any(columns < 1)) {
    stop("inappropriate column parameter. Use integers >=1")
  }
  columns <- round(columns)
  res <- res[, columns]
}




# STEP ONE. Read import and initial processing ---------------------------------------------------------------

code <- readLines(huffmanCode, skipNul = T)


# We split lines into a table (discards the last line, as it apparently doesn't correspond to a byte)
code <- data.table(splitToColumns(code[-length(code)], "\t"))
code[, byte := as.hexmode(as.integer(V3))]


cat("importing reads...\n")

# we ingest all read sequences, ignoring names. This takes memory
reads <- readDNAStringSet(input, use.names = F)

# we place them in a data table
reads <- data.table(SEQ = as.character(reads))

# so we can rapidly discard duplicate read sequences, but retain how many reads
# we had per sequence (much faster than table())
reads <- reads[, .N, by = SEQ]

# forward-sequenced reads should start by AT and end by GC (should be the
# opposite for reverse-sequenced reads)
cat("checking read orientation...\n")
reads[, c("ATstart", "GCend") := .(stri_sub(SEQ, 1, 1) %chin% c("A", "T"), stri_sub(SEQ, 117, 117) %chin% c("G", "C"))]

# contingency table used to evaluate if there are reads to reverse
nums <- matrix(reads[, sum(N), by = .(ATstart, GCend)]$V1, ncol = 2)

# This test should mean there are significantly more apparent reversed reads than expected from sequencing errors
if (chisq.test(nums)$p.value < 0.0001) {
  reads[ATstart == F & GCend == F, SEQ := revCom(SEQ)]
}

reads[, c("ATstart", "GCend") := NULL]




# STAGE TWO. Extraction of read index information -------------------------------------------

cat("extracting read indexing information...\n")

# We extract indexing information in the last 16 bases of a read (excluding the
# very last). We split them into a matrix of individual bases (rows = reads)
IX <- stringToMatrix(stri_sub(reads$SEQ, 101, 116))

IX <- toTrits(IX)

# We get the parity trit by summing the appropriate columns (odd numbers)
parity <- rowSums(IX[, seq(1, 13, 2)]) %% 3L

# We compare it to the 15th trit as an error check. Also, NAs in trits mean errors
error <- parity != IX[, 15] | rowSums(is.na(IX)) > 0
reads <- reads[!error]
IX <- IX[!error, ]

# We obtain file ID from first 2 trits (integer)
fileID <- IX[, 1] * 3L + IX[, 2]

# We obtain read location indices in this file from following 12 trits
idx <- toBase10(IX[, 3:14])


# We add indexing information to the table
reads[, c("FILE", "LOC") := .(fileID, as.integer(idx))]

# we may reclaim some memory
rm(fileID, idx, error, parity, IX)


# to reduce the workload (and free memory), we retain no more than 100 read
# sequences per location index per encoded file. This should be enough to call
# reliable consensuses. We favour sequences without apparent error
cat("sampling best reads...\n")

# values of this column will be TRUE if there are duplicate bases in a read, which shouldn't exist
reads[, err := grepl("AA|CC|GG|TT", SEQ)]

# we sort the table by placing these sequences at the bottom, and favour those represented by many reads
setorder(reads, FILE, LOC, err, -N)

# we retain no more than the first 100 sequences per location per file
reads <- reads[occurrences(stri_c(FILE, LOC, sep = " ")) <= 100, ]



# STAGE TREE, Consensus calling -------------------------------------------------------

# we make a consensus sequences at each location index. We do it before the
# keystream conversion, as a single error may result in a borked read after this
# conversion. We will only convert the (hopefully error-free) consensuses.
cat("generating local consensuses...\n")

# We spot the 100 nt that encode file data into individual bases
reads <- reads[, data.table(stringToMatrix(stri_sub(SEQ, 2, 101)), N, FILE, LOC)]


# The function below counts the number of reads for each of the 4 possible bases at the 100
# positions of a location (takes a base matrix and the number of reads per
# sequences). We will do this for each location index
counts <- function(mat, Ns) {

  # this does it the above for a given type of base
  readCounts <- function(base) {

    # We generate a matrix in which rows represent sequences and columns positions.
    # It contains the number of reads per sequence (identical for all cells of
    # the same row)
    Nmat <- matrix(Ns, nrow(mat), ncol(mat))

    # we remove counts for other bases
    Nmat[mat != base] <- 0L

    # so we can count reads for the base we're interested in
    colSums(Nmat)
  }
  lapply(baseTypes, readCounts)
}


# needed to extract the columns containing the bases from the data table
coln <- names(reads)[1:100]

# We do the base counts by file and location. Note that the ".SD" special symbol
# of data.table would be shorter in code, but it's somehow much slower to run.
# This call creates a file in which each location is represented by 100 rows
# (consensus sequences will appear vertically)
countMat <- reads[, counts(do.call(cbind, mget(coln)), N), by = .(FILE, LOC)]

# We prepare a data.table that will contain consensus sequences
cons <- unique(countMat[, .(FILE, LOC)])

# We turn base counts into a matrix, for operations below
countMat <- as.matrix(countMat[, 3:6, with = F])


# at each position, we get the highest read counts among the 4 bases
maxCounts <- rowMaxs(countMat)

# and corresponding frequency among reads
maxFreqs <- maxCounts / rowSums(countMat)

# we compute a position quality score which takes the base frequency and the
# read count into account. Note that a read count <=1 will results in a score of
# ~0. This may not be appropriate for all experiments.
maxScores <- maxFreqs * (1 - 1 / maxCounts)


# we now select the most frequent base at each position, into a vector
consBases <- baseTypes[rowMatches(maxCounts, countMat)]

# and since there are exactly 100 bases per location, we can easily convert this
# vector to matrix where one row is a location
consBases <- matrix(consBases, ncol = 100, byrow = T)

# we do the same for the scores
consScores <- matrix(maxScores, ncol = 100, byrow = T)
rm(countMat, maxCounts, maxFreqs, maxScores)



# STAGE FOUR, keystream decoding --------------------------------------------------------

# We do the keystream decoding on consensus sequences. The four different keys
# are used as rows of an integer matrix:
cat("keysteam decoding...\n")

ks <- stringToMatrix(c(
  "002000010110102111112122210011122221010102121222022000221201020221002121121000212222021211121122221", 
  "202020122121210120001200210222112020222022222220220001221012111022121120202022211221112202002121022", 
  "221221101200221120220011002222100000020200021121021020122100021201010210202002000101020022121100100", 
  "100122100011112100120210020011102201122122100100120122212000021220022012202201100010212222110222020"))

storage.mode(ks) <- "integer"


# WE convert bases to numbers
d <- chartr("ACGT", "0123", consBases)
storage.mode(d) <- "integer"

# We save the first number (base) of each consensus, for later
first <- d[, 1]

# We compute base-4 differences between adjacent numbers, then subtracting 1
d <- (d[, 2:100] - d[, 1:99]) %% 4L - 1L

# We determine the key to subtract (row of the ks matrix to use) depending on sequence location
phase <- cons$LOC %% 4L + 1L

# We subtract (base 3) the appropriate key and adding 1
d <- (d - ks[phase, ]) %% 3L + 1L

# We sum (base 4) successive trits to obtain base numbers (the reciprocal of the difference step)
d <- cumsumMatrix(cbind(first, d)) %% 4L

# We reverse complement sequences at odd location indices. We use the fact that
# bases are still encoded as numbers, so 3 - base is the complement of base
odd <- cons$LOC %% 2L > 0
d[odd, ] <- 3L - d[odd, ncol(d):1]


# We convert numbers back to bases (chartr() may also work)
d <- matrix(baseTypes[d + 1], nrow = nrow(d))


# We add flattened consensus sequences and average score at each location index to
# the table (we actually don't use these sequences afterwards, but it can be
# used as a visual check)
cons[, c("SEQ", "SCORE") := .(do.call(stri_c, as.data.frame(d)), rowMeans(consScores))]

# we now estimate the length of encoded files, based on the location indices of
# high-quality consensuses. To avoid generating unnecessary long or broken DNA
# sequences

# we do this by creating a table discarding low-score consensuses
temp <- cons[SCORE > 0.75, .(FILE, LOC)]
rows <- 2:nrow(temp)

# We then find "gaps" between successive retained locations (a gap of 6 is arbitrary)
gap <- temp[, c(LOC[rows] - LOC[rows - 1] > 6, F)]

# the last location index with good consensus for each file
ends <- temp[gap == T, min(LOC), by = FILE]$V1

# and the first
starts <- temp[, min(LOC), by = FILE]$V1

# We retain file numbers with enough "good" locations (20 is arbitrary)
valid <- which((ends - starts) > 20) - 1
if (length(valid) == 0) stop("no file with reliable sequence could be found. Exiting.")
rm(temp, starts, gap, rows)
cat(stri_c("file", valid, ", estimated size: ", ends[valid + 1] * 25), sep = "\n")



# STAGE FIVE, global file consensuses generation -----------------------------------------------------

# We now generate the global consensus sequences for files
cat("generating file DNA sequences...\n")

# for this, we cbind the base and score matrices (even if that turns numbers
# into strings). Not very elegant but practical in the following function
d <- cbind(d, consScores)

# We concatenate and overlap successive local consensuses from a file
overlapConsensuses <- function(file) {

  # they will be stored in this matrix, vertically. We have 2*4 columns as we
  # save the base and its score for the 4 overlapping consensuses at each
  # position
  mat <- matrix("", ncol = 8, nrow = ends[file + 1] * 25 + 100)

  # We select the data table rows we need  (that is, excluding what's after the estimated file end)
  f <- cons[, FILE == file & LOC <= ends[file + 1]]

  # and corresponding location indices
  loc <- cons[f == T, LOC]

  # we tolerate that some locations (before the last one) may not be represented, as consensuses overlap
  missing <- setdiff(1:ends[file + 1], loc)
  loc <- c(loc, missing)

  # we extract the relevant rows of the base/score matrix, and add empty rows for missing locations
  bases <- rbind(d[f, ], matrix("", length(missing), 200))

  # phase = the future column of the consensus in the result matrix
  phase <- loc %% 4L

  # this loop concatenates consensuses of the same "phase", since they are adjacent and not overlapping
  for (p in 0:3) {

    # TRUE fo the consensuses of the given phase
    f <- phase == p

    # the first consensus for this phase will start at this offset in the file
    start <- p * 25 + 1

    # the rows of the "d" matrix corresponding to the location indices of this phase
    rows <- which(f)[order(loc[f])]

    # the last position of the concatenated consensus in mat
    end <- start + length(rows) * 100 - 1

    # we can now concatenate consensuses (matrix rows) for each phase and place them in mat
    mat[start:end, p + 1] <- t(bases[rows, 1:100])
    mat[start:end, p + 5] <- t(bases[rows, 101:200])
  }

  # we add file ID as the last column
  cbind(mat, as.character(file))
}


# we apply the function to each valid file. We rbind all in a single matrix
overlaps <- do.call(rbind, lapply(valid, overlapConsensuses))

# we extract the base scores and file ID
baseInfo <- overlaps[, 5:9]
storage.mode(baseInfo) <- "numeric"
baseInfo[is.na(baseInfo)] <- 0


# for each of the 4 bases, we sum its scores over the 4 overlapping consensuses
baseScores <- function(base) {

  # we store scores in a temporary matrix in the function, as we replace them with zeros for the other bases
  scores <- baseInfo[, 1:4]
  scores[overlaps[, 1:4] != base] <- 0
  rowSums(scores)
}


# applying the function for all bases
scores <- vapply(baseTypes, baseScores, numeric(nrow(overlaps)))
maxScores <- rowMaxs(scores)
# hist(maxScores, breaks = 100)

# we retain the base with highest score at each position
consensus <- baseTypes[rowMatches(maxScores, scores)]

# and we create a data table with relevant info (some actually not used afterwards)
consensus <- data.table(
  base = consensus, 
  score = maxScores / rowSums(overlaps[, 1:4] != ""), 
  file = baseInfo[, 5], 
  pos = occurrences(baseInfo[, 5])
  )
# hist(consensus$score, breaks = 100, xlim = 0:1, main = "position quality scores")




# STAGE SIXE, we convert DNA bases to bytes ---------------------------------------------
cat("converting to bytes and saving to file(s)...\n")


# this function does not tolerate missing or erroneous data in files, TO IMPROVE
writeFile <- function(fileID, folder) {

  # the DNA sequence for each file, as a base vector. We prepend base "A", as specified
  bases <- c("A", consensus[file == fileID, base])
  trits <- toTrits(bases)

  # we compute the length of the file, encoded in the first 25 trits.
  len <- toBase10(trits[1:25])
  if (len > length(trits) - 25 | is.na(len)) {
    warning("possible broken file\n")
    len <- length(trits) - 25
  }


  # because "NA" would take 2 chars in the flattened string
  trits[is.na(trits)] <- 3L

  # we flatten the trit vector in a string
  flat <- stri_flatten(trits[26:(len + 25)])

  # we extract 5- and 6-trit words from that string. We need to find where the 6-trit words are
  # all the possible 5-trit words in the string
  word5 <- stri_sub(flat, 1:(len - 4), length = 5)

  # same for 6-trit words
  word6 <- stri_sub(flat, 1:(len - 5), length = 6)

  # we make sure both are in equal numbers (needed to create the matrix below)
  word6 <- c(word6, rep("", length(word5) - length(word6)))

  # we determine whether words are legal, as a 2-column boolean matrix
  allowed <- cbind(word5 %chin% code$V4, word6 %chin% code$V4)

  # to iterate over this matrix (see loop below),
  col <- 1:2
  row <- 1

  # cells of this matrix will be TRUE for the words we retain
  retained <- matrix(F, length(word5), 2)

  # we iterate from row 1
  while (row <= nrow(allowed)) {
    if (allowed[row, col[1]]) {

      # if word is legal we mark its location as TRUE
      retained[row, col[1]] <- T

      # and move to the row of the next word to check (based on word size)
      row <- row + 4 + col[1]
    } else {

      # else if word of alternate length is illegal as well, we abort
      if (!allowed[row, col[2]]) stop(paste("error in file", fileID, "at position", row))

      # if not, we go check the word of alternate length at this position, in the other column
      col <- rev(col)
    }
  }
  
  words <- rbind(word5, word6)[t(retained)]
  bytes <- code[chmatch(words, V4), byte]
  writeBin(as.raw(bytes), con = gsub("//", "/", stri_c(folder, "/file", fileID), fixed = T))
  
  invisible(NULL)
}

for (file in valid) writeFile(file, outputFolder)
