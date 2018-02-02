#I know I use "=" instead of "<-". Don't get mad. 

require(stringi)
require(data.table)	
require(Biostrings)
require(matrixStats)
options("stringsAsFactors"=F)

paths = commandArgs(T)
input = paths[1]					#path to fasta file
outputFolder = paths[2]					#output directory

if(length(paths) < 3) {
	huffmanCode = gsub("//","/", stri_c(outputFolder, "/View_huff3.cd.new.correct"), fixed =T)
} else {
	huffmanCode = paths[3]
}

paths = c(input, outputFolder, huffmanCode)
missing = !file.exists(paths)

if(any(missing)) stop(paste("missing file(s):", paste(paths[missing], collapse =", ")))


###########defining "basic" functions" used in the script

baseTypes = c("A","C","G","T")			

conv = matrix(c(NA, 0:2, 2L, NA, 0:1, 1:2, NA, 0L, 0:2, NA), 4, byrow = T)	#matrix used to convert successive bases to trits, used below
colnames(conv) = baseTypes; rownames(conv) = baseTypes

toTrits = function(mat) {							#converts a matrix of DNA bases (sequences as rows) to integer trits according to Goldman et al's scheme
	if(!is.matrix(mat)) mat = rbind(mat)					#in case mat is a vector, we convert it to a 1-row matrix (otherwise we would have dimension issues)
	vapply(2:ncol(mat), function(col) conv[cbind(mat[,col-1], mat[,col])], integer(nrow(mat)))
}

stringToMatrix = function(vector) {						#splits a string vector into a matrix with one character per cell. Words are rows. The number of columns is the length of the longest word. We use as substring function (stri_sub()), which appears to be faster than a split
	nchars = stri_length(vector)
	m = max(nchars)
	mat = matrix(NA,nrow = length(vector), ncol = m,byrow=T)
	for (i in 1:m) {
		mat[,i] = stri_sub(vector,i,i)		
	}
	mat
}

revCom = function(seqs) {  							#reverse complements a character vector representing DNA sequences. Ambiguities are not treated
	stri_reverse(chartr("acgtACGT","tgcaTGCA", seqs))
}


toBase10 = function(mat, base=3) {						#converts numbers of a given base to base 10, assuming digits are in separate columns of a matrix, so that rows represent numbers read from left to right
	if(!is.matrix(mat)) mat = rbind(mat)	
	res = 0									
	ncol = ncol(mat)
	for (col in ncol:1) {
		res = res + mat[,col]*base^(ncol-col)
	}
	res
}

occurences = function(vect) {							#tells the occurence of the element found, for all positions of a vector (e.g. 1 if it's the first time it appears, 2 if the 2nd time, etc)
	dt = data.table(vect, pos = 1:length(vect))
	counts = dt[,.(occ=1:.N, pos), by = vect]
	res = integer(length(vect))
	res[counts$pos] = counts$occ
	res
}

rowMatches = function(values, mat) {						#returns the first position (column) of each element of "values" in each corresponding row of matrix "mat" (recycling may work, as values can have 1 element)
	res = which(mat==values, arr.ind = T)
	res[match(1:nrow(mat), res[,1]), 2]
}

cumsumMatrix = function(mat) {							#returns a matrix in which colums represent the cumulated sums of the colums of mat, starting from the left
	mat = cbind(mat, -rowSums(mat))
	cs = matrix(cumsum(t(mat)), nrow(mat), byrow = T)
	cs[,-ncol(cs)]
}

splitToColumns = function(vector, split, columns) {  				#splits a character vector into columns (like excel). A column range can be specified as an integer vector and "cells" are filled with NAs if needed. By default, it generate as many columns as necessary given the character(s) used to split the vector elements.
	vector = as.character(vector)
	res = stri_split_fixed(vector,split, omit_empty = NA, simplify = T)
	if (missing(columns)) {
		columns = 1:ncol(res)
	}
	if (any(columns < 1)) {
		stop("inappropriate column parameter. Use integers >=1")
	}
	columns =round(columns)
	res = res[,columns]
}

###################################################
	       
code = readLines(huffmanCode, skipNul = T)					#imports the huffman code. read.table() doesn't like it.

code = data.table(splitToColumns(code[-length(code)] ,"\t"))			#split lines into a table (discards the last line, as it apparently doesn't correspond to a byte)
code[,byte := as.hexmode(as.integer(V3))]					
	       
cat("importing reads...\n")
reads = readDNAStringSet(input, use.names = F)  				#we ingest all read sequences, ignoring names. This takes memory
reads = data.table(SEQ = as.character(reads))					#we place them in a data table
reads = reads[,.N, by = SEQ]							#so we can rapidly discard duplicate read sequences, but retain how many reads we had per sequence (much faster than table())

#forward-sequenced reads should start by AT and end by GC (should be the opposite for reverse-sequenced reads)
cat("checking read orientation...\n")
reads[,c("ATstart","GCend") := .(stri_sub(SEQ, 1, 1) %chin% c("A","T"), stri_sub(SEQ, 117, 117) %chin% c("G","C"))]	
nums = matrix(reads[, sum(N), by = .(ATstart, GCend)]$V1, ncol = 2)		#contigency table used to evaluate if there are reads to reverse 
if (chisq.test(nums)$p.value < 0.0001) {											#should mean there are significantly more apparent reversed reads than expected from sequencing errors
	reads[ATstart == F & GCend == F, SEQ := revCom(SEQ)]
}

reads[,c("ATstart", "GCend") := NULL]

cat("extracting read indexing information...\n")
IX = stringToMatrix(stri_sub(reads$SEQ, 101, 116))				#extracts indexing information in the last 16 bases of a read (excluding the very last). We split them into a matrix of individual bases (rows = reads)

IX = toTrits(IX)												
parity = rowSums(IX[,seq(1, 13, 2)]) %% 3L 					#gets the parity trit by summing the appropriate columns (odd numbers)
error = parity != IX[,15] | rowSums(is.na(IX)) > 0				#compares it to the 15th trit as an error check. Also, NAs in trits mean errors
reads = reads[!error]; IX = IX[!error,]						
fileID = IX[,1]*3L + IX[, 2]							#obtains file ID from first 2 trits (integer)
idx= toBase10(IX[,3:14])							#obtains read location indices in this file from following 12 trits

reads[,c("FILE", "LOC") := .(fileID, as.integer(idx))]				#adds indexing information to the table 
rm(fileID, idx, error, parity, IX)						#we may reclame some memory

#to reduce the workload (and free memory), we retain no more than 100 read sequences per location index per encoded file. This should be enough to call reliable consensuses. We favor sequences without apparent error
cat("sampling best reads...\n")
reads[,err := grepl("AA|CC|GG|TT", SEQ)]					#values of this column will be TRUE if there are duplicate bases in a read, which shouldn't exist
setorder(reads, FILE, LOC, err, -N)						#we sort the table by placing these sequences at the bottom, and favor those represented by many reads
reads = reads[occurences(stri_c(FILE, LOC, sep = " ")) <= 100, ]		#we retain no more than the first 100 sequences per location per file

#now making a consensus sequences at each location index. We do it before the keystream conversion, as a single error may result in a borked read after this conversion. We will only convert the (hopefully error-free) consensuses.
cat("generating local consensuses...\n")
reads = reads[, data.table(stringToMatrix(stri_sub(SEQ, 2, 101)), N, FILE, LOC)]	#splits the 100 nt that encode file data into individual bases

counts = function(mat, Ns) {								#this counts the number of reads for each of the 4 possible bases at the 100 positions of a location (takes a base matrix and the number of reads per sequences). We will do this for each location index
	readCounts = function(base) {							#this does it the above for a given type of base
		Nmat = matrix(Ns, nrow(mat), ncol(mat))					#generate a matrix in which rows represent sequences and columns positions. It contains the number of reads per sequence (identical for all cells of the same row)
		Nmat[mat != base] = 0L							#we remove counts for other bases
		colSums(Nmat)								#so we can count reads for the base we're interested in
	}
	lapply(baseTypes, readCounts)
}

coln = names(reads)[1:100]								#needed to extract the columns containing the bases from the data table
countMat = reads[,counts(do.call(cbind, mget(coln)), N), by = .(FILE, LOC)]		#doing the base counts by file and location. Note that the ".SD" special symbol of data.table would be shorter in code, but it's somehow much slower to run. This call creates a file in which each location is represented by 100 rows (consensus sequences will appear vertically)
cons = unique(countMat[, .(FILE, LOC)])							#preparing a data.table that will contain consensus sequences
countMat = as.matrix(countMat[,3:6, with = F])						#turning base counts into a matrix, for operations below 

maxCounts = rowMaxs(countMat)								#at each position, we get the highest read counts among the 4 bases
maxFreqs = maxCounts/rowSums(countMat)							#and corresponding frequency among reads
maxScores = maxFreqs *(1 - 1/maxCounts)							#we compute a position quality score which takes the base frequency and the read count into account. Note that a read count <=1 will results in a score of ~0. This may not be appropriate for all experiments. 

consBases = baseTypes[rowMatches(maxCounts, countMat)]					#we now select the most frequent base at each position, into a vector
consBases = matrix(consBases, ncol = 100, byrow = T)					#and since there are exactly 100 bases per location, we can easily convert this vector to matrix where one row is a location
consScores = matrix(maxScores, ncol = 100, byrow = T)					#we do the same for the scores
rm(countMat, maxCounts, maxFreqs, maxScores)

#doing the keystream decoding on consensus sequences. The four different keys are used as rows of an integer matrix:
cat("keystream decoding...\n")
ks = stringToMatrix(c("002000010110102111112122210011122221010102121222022000221201020221002121121000212222021211121122221", "202020122121210120001200210222112020222022222220220001221012111022121120202022211221112202002121022", "221221101200221120220011002222100000020200021121021020122100021201010210202002000101020022121100100", "100122100011112100120210020011102201122122100100120122212000021220022012202201100010212222110222020"))
storage.mode(ks) = "integer"

d = chartr("ACGT","0123", consBases)							#converts bases to numbers
storage.mode(d) = "integer"
first = d[,1]										#saves the first number (base) of each consensus, for later
d = (d[,2:100] - d[,1:99]) %% 4L - 1L							#computing base-4 differences between adjacent numbers, then substracting 1
phase = cons$LOC %% 4L + 1L								#determines the key to substract (row of the ks matrix to use) depending on sequence location
d = (d-ks[phase,]) %% 3L + 1L								#substracting (base 3) the appropriate key and adding 1
d = cumsumMatrix(cbind(first, d)) %% 4L							#summing (base 4) successive trits to obtain base numbers (the reciprocal of the difference step)

#reverse complements sequences at odd location indices. We use the fact that bases are still encoded as numbers, so 3 - base is the complement of base
odd = cons$LOC %% 2L > 0
d[odd, ] = 3L - d[odd, ncol(d):1]

d = matrix(baseTypes[d+1], nrow=nrow(d))						#converting numbers back to bases (chartr() may also work)

cons[, c("SEQ", "SCORE") := .(do.call(stri_c, as.data.frame(d)), rowMeans(consScores))]			#adds flattened consensus sequences and average score at each location index to the table (we actually don't use these sequences afterwards, but it can be used as a visual check)

#we now estimate the length of encoded files, based on the location indices of high-quality consensuses. To avoid generating unnecessary long or broken DNA sequences
temp = cons[SCORE > 0.75, .(FILE, LOC)]							#we do this by creating a table discarding low-score consensuses
rows = 2:nrow(temp)
gap = temp[, c(LOC[rows] - LOC[rows-1] > 6, F)]						#then finding "gaps" between successive retained locations (a gap of 6 is arbitrary)
ends = temp[gap == T, min(LOC), by = FILE]$V1						#the last location index with good consensus for each file
starts = temp[, min(LOC), by = FILE]$V1							#and the first
valid = which((ends - starts) > 20) - 1							#retains file numbers with enough "good" locations (20 is arbitrary)
if(length(valid)==0) stop("no file with reliable sequence could be found. Exiting.")
rm(temp, starts, gap, rows)
cat(stri_c("file", valid, ", estimated size: ", ends[valid+1]*25), sep="\n")


#now generating the global consensus sequences for files
cat("generating file DNA sequences...\n")
d = cbind(d, consScores)								#for this, we cbind the base and score matrices (even if that turns numbers into strings). Not very elegant but practical in the following function
overlapConsensuses = function(file) {							#concatenates and overlaps successive local consensuses from a file
	mat = matrix("", ncol= 8, nrow = ends[file+1]*25+100)				#they will be stored in this matrix, vertically. We have 2*4 columns as we save the base and its score for the 4 overlapping consensuses at each position 
	f = cons[,FILE == file & LOC <= ends[file+1]]					#selects the data table rows we need  (that is, excluding what's after the estimated file end)
	loc = cons[f==T, LOC]								#and corresponding location indices	
	missing = setdiff(1:ends[file+1], loc)						#we tolerate that some locations (before the last one) may not be represented, as consensuses overlap
	loc = c(loc, missing)									
	bases = rbind(d[f,], matrix("", length(missing), 200))				#we extract the relevant rows of the base/score matrix, and add empty rows for missing locations
	phase = loc %% 4L								#phase = the future column of the consensus in the result matrix
	for (p in 0:3) {								#this loop concatenates consensuses of the same "phase", since they are adjacent and not overlapping
		f = phase == p								#TRUE fo the consensuses of the given phase
		start = p*25+1								#the first consensus for this phase will start at this offset in the file
		rows = which(f)[order(loc[f])]						#the rows of the "d" matrix corresponding to the location indices of this phase
		end = start + length(rows)*100-1					#the last position of the concatenated consensus in mat
		mat[start:end, p+1] = t(bases[rows,1:100])				#we can now concatenate consensuses (matrix rows) for each phase and place them in mat
		mat[start:end, p+5] = t(bases[rows,101:200])
	}	
	cbind(mat, as.character(file))							#adds file ID as the last column
}

overlaps = do.call(rbind, lapply(valid, overlapConsensuses))				#applying the function to each valid file. We rbind all in a single matrix
baseInfo = overlaps[,5:9]; storage.mode(baseInfo) = "numeric"				#extracts the base scores and file ID
baseInfo[is.na(baseInfo)] = 0

baseScores = function(base) {								#for each of the 4 bases, we sum its scores over the 4 overlapping consensuses
	scores = baseInfo[,1:4]								#we store scores in a temporary matrix in the function, as we replace them with zeros for the other bases
	scores[overlaps[,1:4] != base] = 0
	rowSums(scores)
}

scores = vapply(baseTypes, baseScores, numeric(nrow(overlaps)))				#applying the function for all bases
maxScores = rowMaxs(scores)
#hist(maxScores, breaks = 100)
consensus = baseTypes[rowMatches(maxScores, scores)]					#we retain the base with highest score at each position
consensus = data.table(base = consensus, score = maxScores/rowSums(overlaps[,1:4] != ""), file = baseInfo[,5], pos = occurences(baseInfo[,5]))	#and create a data table with relevant info (some actually not used afterwards)
#hist(consensus$score, breaks = 100, xlim = 0:1, main = "position quality scores")

#now converting DNA to bytes
cat("converting to bytes and saving to file(s)...\n")

writeFile = function(fileID, folder) {							#this function does not tolerate missing or erronenous data in files, TO IMPROVE
	bases = c("A", consensus[file ==  fileID, base])				#the DNA sequence for each file, as a base vector. We prepend base "A", as specified
	trits = toTrits(bases)
	len = toBase10(trits[1:25])							#computes the length of the file, encoded in the first 25 trits.
	if (len > length(trits)-25 | is.na(len)) {
			warning("possible broken file\n")
			len = length(trits)-25
	}
	
	trits[is.na(trits)] = 3L							#because "NA" would take 2 chars in the flattened string
	flat = stri_flatten(trits[26:(len+25)])						#flattens the trit vector in a string
	
	#extracting 5- and 6-trit words from that string. We need to find where the 6-trit words are
	word5 = stri_sub(flat, 1:(len-4), length = 5)					#all the possible 5-trit words in the string
	word6 = stri_sub(flat, 1:(len-5), length = 6)					#same for 6-trit words
	word6 = c(word6, rep("", length(word5)-length(word6)))				#makes sure both are in equal numbers (needed to create the matrix below)
	allowed = cbind(word5 %chin% code$V4, word6 %chin% code$V4)			#tells whether words are legal, as a 2-column boolean matrix
	col = 1:2 ; row = 1 								#to iterate over this matrix (see loop below), 
	retained = matrix(F,length(word5), 2)						#cells of this matrix will be TRUE for the words we retain
	while (row <= nrow(allowed)) {							#iterates from row 1
		if(allowed[row,col[1]]) {								
			retained[row, col[1]] = T					# if word is legal we mark its location as TRUE
			row = row + 4 + col[1]						#and move to the row of the next word to check (based on word size)
		} else {												
			if(!allowed[row, col[2]]) stop(paste("error in file", fileID, "at position", row))		#else if word of alternate length is illegal as well, we abort 
			col = rev(col)							#if not, we go check the word of alternate length at this position, in the other column
		}
	}
	words = rbind(word5, word6)[t(retained)]	
	bytes = code[chmatch(words, V4), byte]
	writeBin(as.raw(bytes), con = gsub("//","/", stri_c(folder, "/file", fileID), fixed = T))
	invisible(NULL)
}

for (file in valid) writeFile(file, outputFolder)

