# seqDataTools.R


detectFastqReadFormat <- function( filein) {

	# get ready to read the fastq file in chuncks...
	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")

	# there is often garbage at the front...
	# so step past it...
	tmp <- readLines( con=conIn, n=10000)

	# now use the last few
	nlines <- length(tmp)
	tmp <- tmp[ max( c(nlines-40+1, 1)):nlines]
	nlines <- length(tmp)
	close( conIn)

	readids <- tmp[ seq( 1, nlines, by=4)]
	readseqs <- tmp[ seq( 2, nlines, by=4)]
	scores <- tmp[ seq( 4, nlines, by=4)]

	# verify bases...
	mybases <- base::unlist( strsplit( readseqs, split=""))
	whobad <- which( !( mybases %in% c("A","C","G","T","N",".")))
	if ( length(whobad) > 0) {
		stop( paste( "Bad fastq file:  non-DNA bases detected:\nFile: ", filein, 
				"\nBad Bases: ", mybases[whobad]))
	}

	# see what separators are used in the read IDs
	readIDtype <- detectReadIDformat( readids)

	# detect the score encryption
	try33 <- phredScoreStringToInt( scores, "Phred33")
	try64 <- phredScoreStringToInt( scores, "Phred64")
	# none should be negative
	scoreType <- "Phred64"
	nTooSmall33 <- sum( try33 < 1)
	nTooSmall64 <- sum( try64 < 1)
	if ( nTooSmall33 < nTooSmall64) {
		scoreType <- "Phred33"
	}

	out <- list( "readIDtype"=readIDtype, "scoreType"=scoreType)
	return( out)
}


`detectReadIDformat` <- function( readids) {

	N <- length( readids)

	# see what separators are used in the read IDs
	colons <- base::unlist( gregexpr( ":", readids, fixed=T))
	underscores <- base::unlist( gregexpr( "_", readids, fixed=T))
	periods <- base::unlist( gregexpr( ".", readids, fixed=T))
	idealSolexa1.3 <- length(colons) / (4 * N)
	idealRosetta <- length(underscores) / (3 * N)
	idealSolexa1.8 <- length(colons) / (9 * N)
	idealUnusual1 <- length(periods) / (5 * N)
	
	bestGuess <- c( idealSolexa1.3, idealRosetta, idealSolexa1.8, idealUnusual1)
	names(bestGuess) <- c( "Solexa_1.3", "Rosetta", "Solexa_1.8", "Unusual1")
	best <- which.min( abs( bestGuess - 1))
	readIDtype <- names(bestGuess)[best]

	return( readIDtype)
}


extractReadIDlaneNumber <- function( txt, isFastq=TRUE, readIDtype="Rosetta") {

	# the read ID is " @<lane>_<tile>_<x>_<y> "
	out <- rep( 0, times=length( txt))

	if ( readIDtype == "Rosetta") {
		seps <- gregexpr( "_", txt[1], fixed=TRUE)[[1]]
		nsep <- length(seps) 
		sep <- seps[ nsep - 2]
		out <- as.integer( base::substr( txt, (sep-1), (sep-1)))
	} else if ( readIDtype == "Solexa_1.3") {
		# in derisi IDs there is a ":n:m:x:y" field...
		seps <- gregexpr( ":", txt[1], fixed=TRUE)[[1]]
		sep <- if (length(seps) > 0) seps[1] else 0
		if (length(seps) > 4) sep <- seps[ length(seps) - 3]
		out <- ifelse( sep > 0, as.integer( base::substr( txt, (sep+1), (sep+1))), 0)
	} else if ( readIDtype == "Solexa_1.8") {
		# in newest IDs there is a "instrumentID:runID:flowCellID:n:m:x:y a:b:c:" field...
		# in derisi IDs there is a ":n:m:x:y" field...
		seps <- gregexpr( ":", txt[1], fixed=TRUE)[[1]]
		out <- ifelse( length(seps) > 3, as.integer( base::substr( txt, (seps[3]+1), (seps[4]-1))), 0)
	} else if ( readIDtype == "Unusual1") {
		# in IDs there is a ".n.m.x.y.readnumber" field...
		seps <- gregexpr( ".", txt[1], fixed=TRUE)[[1]]
		out <- ifelse( length(seps) > 1, as.integer( base::substr( txt, (seps[1]+1), (seps[2]-1))), 0)
	} else {
		stop( paste( "extractReadIDlaneNumber:  unknowm readID type: ", readIDtype, 
				"Known:  Rosetta, Solexa_1.3, Solexa_1.8"))
	}

	if ( ! all( out > 0)) {
		cat("\nExtracting ReadID Lane Number failed:   readID sample=\n")
		print( head( txt))
	}

	return( out)
}


extractReadIDterms <- function( txt, isFastq=TRUE, readIDtype="Rosetta") {

	# the read ID is " @<lane>_<tile>_<x>_<y> " for older Rosetta style
	# in DeRisi and newer Solexa, the  ID is a "@name:lane:tile:x:y" style
	# in newest Solexa, the  ID is a "@name:str::lane:tile:x:y a:b:c:" style
	outLane <- outTile <- outX <- outY <- rep( 0, times=length( txt))

	strIn <- txt

	skipChar <- 0
	if ( readIDtype == "Rosetta") {
		if ( isFastq) skipChar <- 1
		splitMark <- "_"
		wantedTerms <- 1:4
		fixed <- TRUE
	} else if ( readIDtype == "Solexa_1.3") {
		splitMark <- ":"
		wantedTerms <- 2:5
		fixed <- TRUE
	} else if ( readIDtype == "Solexa_1.8") {
		splitMark <- ":| "
		wantedTerms <- 4:11
		fixed <- FALSE
	} else if ( readIDtype == "Unusual1") {
		splitMark <- "."
		wantedTerms <- 2:5
		fixed <- TRUE
	} else {
		stop( paste( "extractReadIDterms:  unknowm readID type: ", readIDtype, 
				"Known:  Rosetta, Solexa_1.3, Solexa_1.8"))
	}
	markPtrs <- gregexpr( splitMark, strIn[1], fixed=fixed)[[1]]
	nTermsEach <- length( markPtrs) + 1

	# solexa_1.8 does not always have a barcode
	hasBarcode <- FALSE
	if ( readIDtype == "Solexa_1.8") {
		hasBarcode <- TRUE
		if ( markPtrs[ length(markPtrs)] == base::nchar( strIn[1])) {
			nTermsEach <- nTermsEach - 1
			hasBarcode <- FALSE
		}
		# some datasets will have some pieces trimmed away already
		if( max( wantedTerms) > nTermsEach) {
			wantedTerms <- 4:nTermsEach
			hasBarcode <- FALSE
		}
	}
	# some 1.3 readIDs have colons in the prefix before the lane,tile,etc. part
	if ( readIDtype == "Solexa_1.3") {
		if ( length(markPtrs) > 4) {
			nn <- length(markPtrs)
			wantedTerms <- (nn-4+2) : (nn+1)
		}
	}
	# some Rosetta readIDs have text before the lane number
	if ( readIDtype == "Rosetta") {
		if ( length(markPtrs) > 3) {
			nn <- length(markPtrs)
			wantedTerms <- (nn-2) : (nn+1)
		}
	}


	allTerms <- strsplit( strIn, split=splitMark, fixed=fixed)
	mTerms <- matrix( base::unlist(allTerms), nrow=nTermsEach, ncol=length(strIn))
	if (skipChar > 0 && wantedTerms[1] == 1) {
		mTerms[1, ] <- base::substr( mTerms[1, ], skipChar+1, 10)
	}
	outLane <- as.integer( mTerms[ wantedTerms[1], ])
	if ( any( is.na(outLane))) {
		cat("\nFound non-integer Lane Numbers.    Tried ID field: ", wantedTerms[1])
		cat( "\nEncounterd:  ",  mTerms[ wantedTerms[1], head( which( is.na(outLane)))])
	}
	outTile <- as.integer( mTerms[ wantedTerms[2], ])
	if ( any( is.na(outTile))) {
		cat("\nFound non-integer Tile Numbers.    Tried ID field: ", wantedTerms[2])
		cat( "\nEncounterd:  ",  mTerms[ wantedTerms[2], head( which( is.na(outTile)))])
	}
	outX <- as.integer( mTerms[ wantedTerms[3], ])
	if ( any( is.na(outX))) {
		cat("\nFound non-integer X_position Values.    Tried ID field: ", wantedTerms[3])
		cat( "\nEncounterd:  ",  mTerms[ wantedTerms[3], head( which( is.na(outX)))])
	}
	# Y term may have junk...
	if ( readIDtype %in% c( "Rosetta", "Solexa_1.3", "Solexa_1.8")) {
		tmp <- mTerms[ wantedTerms[4], ] 
		tmp <- sub( "( |#|/).*$", "", tmp) 
		outY <- as.integer( tmp)
	} else {
		outY <- as.integer( mTerms[ wantedTerms[4], ])
	}
	if ( any( is.na(outY))) {
		cat("\nFound non-integer Y_position Values.    Tried ID field: ", wantedTerms[4])
		cat( "\nEncounterd:  ",  mTerms[ wantedTerms[4], head( which( is.na(outY)))])
	}

	isFiltered <- barcode <- NULL
	if ( readIDtype %in% c( "Solexa_1.8")) {
		if ( length( wantedTerms) >= 6) isFiltered <- mTerms[ wantedTerms[6], ]
		if ( length( wantedTerms) >= 8 && hasBarcode) barcode <- mTerms[ wantedTerms[8], ]
	}

	if ( any( is.na(outLane))) cat("\nBad Lane Numbers: ", mTerms[ wantedTerms[1], 
			head( which( is.na(outLane)))])
	
	if ( ! all( outLane > 0)) {
		cat("\nExtracting ReadID Terms failed:   readID sample=\n")
		print( head( txt))
	}

	return( list( "Lane"=outLane, "Tile"=outTile, "X"=outX, "Y"=outY, "Filter"=isFiltered,
			"Barcode"=barcode))
}



partitionFastqByReadIDfield <- function( filein, field=c( "Lane", "Filter")) {

	# partition one fastq file of multiple lanes of data into separate .fastq files for each lane.

	field <- match.arg( field)

	fileToUse <- allowCompressedFileName( filein)
	if ( ! file.exists( fileToUse)) stop( paste("Can't find input file: ", fileToUse))
	compressOutput <- ( regexpr( ".gz$", fileToUse) > 0)

	# get the type of readIDs
	ans <- detectFastqReadFormat( fileToUse)
	readIDtype <- ans$readIDtype
	cat( "\nFastq ReadID format:  ", readIDtype, "\n")

	conIn <- openCompressedFile( fileToUse, open="r")

	chunkSize <- 400000
	nread <- 0
	fieldSet <- countSet <- fileSet <- vector()
	conSet <- list()
	partitionOverhead <- list( "fields"=fieldSet, "cons"=conSet, "counts"=countSet, "files"=fileSet)

	repeat {
		chunk <- readLines( conIn, n=chunkSize)
		if ( length( chunk) < 1) break

		nread <- nread + length( chunk)
		cat( "\rN_Reads: ", formatC( as.integer(nread/4), big.mark=","))
		partitionOverhead <- partitionFastq( chunk, field=field, readIDtype=readIDtype, 
				overhead=partitionOverhead, basefile=basename( fileToUse), 
				compressOutput=compressOutput)
		cat ( "  ", partitionOverhead$fields, "   ", formatC( as.integer( 
				partitionOverhead$counts), big.mark=","))
	}

	close( conIn)
	cat( "\nN_reads in original:  ", formatC( as.integer(nread/4), big.mark=","))
	for( i in 1:length( partitionOverhead$fields)) {
		thisField <- partitionOverhead$fields[ i]
		con <- partitionOverhead$cons[[ i]]
		cnt <- partitionOverhead$counts[ i]
		fil <- partitionOverhead$files[ i]
		close( con)
		cat( "\n",field, "=", thisField, "\tWrote file: ", fil, "\tN_Reads: ", formatC( as.integer(cnt),
				big.mark=","))
	}
	cat( "\n")
	return()
}


partitionFastq <- function( txt, field="Lane", readIDtype="Rosetta", overhead, basefile,
			compressOutput=FALSE) {

	# given a chunk of .fastq file text, use a field in the id to partition...

	# the id line is 1st, and it has the "lane" number in it...
	idsTxt <- txt[ seq( 1, length(txt), by=4)]
	idTerms <- extractReadIDterms( idsTxt, readIDtype=readIDtype)

	if ( ! (field %in% names(idTerms))) {
		cat( "\nInvalid ReadID field:  ", field, 
			"\nChoices are:   ", names( idTerms), "\n")
		stop()
	}

	# factor by lane number
	myFields <- idTerms[[ field]]
	fieldFac <- factor( myFields)
	thesefields <- levels( fieldFac)

	# any new lanes get their files, etc., set up...
	for ( i in 1:length(thesefields)) {
		thisfield <- thesefields[i]
		if ( thisfield %in% overhead$fields) next
		# this is a new lane, add it, open a new connection, etc.
		nnow <- length( overhead$fields) + 1
		overhead$fields[nnow] <- thisfield

		if ( regexpr( ".fastq", basefile, fixed=TRUE) > 0) {
			fileout <- sub( ".fastq", paste( ".", field, thisfield, ".fastq", sep=""), basefile, fixed=TRUE)
		} else if ( regexpr( ".fq", basefile, fixed=TRUE) > 0) {
			fileout <- sub( ".fq", paste( ".", field, thisfield, ".fq", sep=""), basefile, fixed=TRUE)
		} else {
			fileout <- paste( basefile, ".", field, thisfield, sep="")
		}
		if (compressOutput) {
			if ( regexpr( ".gz$", fileout) < 1) {
				fileout <- paste( fileout, "gz", sep=".")
			}
			con <- gzfile(fileout, open="wb")
		} else {
			con <- file(fileout, open="wt")
		}
		overhead$files[nnow] <- fileout
		overhead$cons[[nnow]] <- con
		overhead$counts[nnow] <- 0
	}

	# now partition those lines of text
	for ( field in thesefields) {
		whichfield <- base::match( field, overhead$fields)
		mycon <- overhead$cons[[ whichfield]]
		mycnt <- overhead$counts[ whichfield]

		who <- which( myFields == field)
		# turn these read locations into the actual 4 line pointers for each read
		starts <- ((who-1) * 4) + 1
		ends <- starts + 3
		lines <- base::unlist( base::mapply( FUN=`:`, starts, ends))

		# write those lines to that file
		writeLines( txt[ lines], con=mycon)

		# update the overhead
		overhead$counts[ whichfield] <- mycnt + length( starts)
	}

	return( overhead)
}


seq2fastq <- function( filein, fileout, phred=TRUE, split="\t") {

	# turn a .seq file from Rosetta ( compressed or not) having Solexa-type quality scores
	# into a .fastq file with Phred-type quality scores ( compressed or not)

	filein <- allowCompressedFileName(filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")

	if ( regexpr( ".gz$", fileout) > 0) {
		conOut <- gzfile( fileout, open="w")
	} else {
		conOut <- file( fileout, open="w")
	}

	cat( "\nInput File: ", filein)

	chunkSize <- 100000
	nread <- 0

	repeat {
		chunk <- readLines( conIn, n=chunkSize)
		if ( length( chunk) < 1) break

		nread <- nread + length( chunk)
		newchunk <- illumina2Fastq( chunk, phred=phred, split=split)
		writeLines( newchunk, conOut)
		cat( "  ",formatC(nread, format="d", big.mark=","))
	}

	close( conIn)
	close( conOut)
	cat( "\nN_reads converted: ", nread, "\n")
}


illumina2Fastq <- function( txt, phred=TRUE, split=":") {

	# given lines of text from a illumina .SEQ file of raw reads,
	# transform it into a .FASTQ file for input to an alignment program
	terms <- strsplit( txt, split=split, fixed=TRUE)

	# we should trap any bad records before we 'matrix' them...
	nExpect <- 7
	if ( split == "\t") nExpect <- 11
	nGot <- sapply( terms, length)
	if ( (nUse <- median( nGot)) > nExpect) nExpect <- nUse

	badOnes <- which( nGot < nExpect)
	if ( length( badOnes) > 0) {
		terms <- terms[ -badOnes]
		badtxt <- txt[ badOnes]
		cat( "\nCorrupt Lines:   skipping ", length(badOnes),"\n")
		print( badtxt)
	}

	flat <- base::unlist( terms)
	termsM <- matrix( flat, nrow=nExpect, ncol=length( terms))

	# older Illumina style was 'colon' separated...
	if ( split == ":") {
	
	# first field is id, ignore for now
	# next four are "lane / tile / X / Y"
	myID <- base::paste( termsM[2,], termsM[3,], termsM[4,], termsM[5,], sep="_")
	laneID <- termsM[ 2, ]

	# next is the RNA read string
	mySEQ <- termsM[ 6, ]

	# last is the quality scores, as space-separated "Solexa-type" integer scores
	# OR as a "Solexa-type" character string...
	myScores <- termsM[ 7, ]

	# try to detect which...
	hasSpaces <- regexpr( " ", myScores, fixed=TRUE)
	if ( all( hasSpaces > 0)) {
		# nothing to do....
	} else {
		fac <- factor( myScores)
		uniqueScores <- levels( fac)
		facPtrs <- tapply( 1:length(myScores), fac, FUN=NULL)
		cat( "  Solexa Conv Speedup=", formatC( (length(myScores)/length(uniqueScores)), digits=2, format="f"))
		uniqueSolexaStrs <- sapply( uniqueScores, FUN=solexaCharScoreToSolexaIntScore)
		myScores <- uniqueSolexaStrs[ facPtrs]
	}
	# convert to Phred style scores...
	if ( phred) {
		fac <- factor( myScores)
		uniqueScores <- levels( fac)
		facPtrs <- tapply( 1:length(myScores), fac, FUN=NULL)
		cat( "  Phred Conv Speedup=", formatC( (length(myScores)/length(uniqueScores)), digits=2, format="f"))
		uniquePhredStrs <- sapply( uniqueScores, FUN=solexaToPhred)
		myScores <- uniquePhredStrs[ facPtrs]
	}
	}

	if ( split == "\t") {  #if split is by tab, newer format with different spacing, etc...

	# next four are "lane / tile / X / Y"
	myIDa <- base::paste( termsM[1], termsM[2,], sep="_") 
	myIDb <- base::paste( termsM[3,], termsM[4,], termsM[5,], termsM[6,], sep=":")
	myIDc <- base::paste( termsM[7], termsM[8,], sep="/") 
	myID <- base::paste( myIDa, ":", myIDb, "#", myIDc, sep="")
	laneID <- termsM[ 3, ]

	# next is the RNA read string
	mySEQ <- termsM[ 9, ]

	# "Solexa-type" character string...
	myScores <- termsM[ 10, ]
	}


	# make the fastq format lines...
	out <- base::paste( "@", myID, "\n", mySEQ, "\n+\n", myScores, sep="")
	return( out)
}


fastqPatternSearchDetails <- function( filein, pattern="AGAGATCGGAAGATCTCGTATGCC") {

	# find all reads that contain a pattern string

	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))

	ans <- detectFastqReadFormat( filein)
	readIDtype <- ans$readIDtype

	conIn <- openCompressedFile( filein, open="r")

	chunkSize <- 400000
	nreads <- nhits <- 0

	bigTileYes <- bigXYes <- bigYYes <- vector()
	bigTileNo <- bigXNo <- bigYNo <- vector()

	repeat {
		txt <- readLines( conIn, n=chunkSize)
		if ( length( txt) < 1) break

		nread <- round(length( txt)/4)
		nreads <- nreads + nread
		cat( "\nN_Reads: ", prettyNum( as.integer(nreads), big.mark=","))

		# extract all the lane,tile,X,Y info from the read IDs
		idsTxt <- txt[ seq( 1, length(txt), by=4)]
		idTerms <- extractReadIDterms( idsTxt, readIDtype=readIDtype)

		# turn the reads into a searchable set, and count how often we see the pattern
		readTxt <- txt[ seq( 2, length(txt), by=4)]
		subject <- DNAStringSet( readTxt)

		cnts <- vcountPattern( pattern, subject, max.mismatch=1)
		hits <- which( cnts > 0)
		nhit <- length( hits)
		nhits <- nhits + nhit
		cat( "  N_hits: ", nhit, "  Percent Hits: ", as.percent( nhit, big.value=nread))

		bigTileYes <- append( bigTileYes, idTerms$Tile[ hits])
		bigXYes <- append( bigXYes, idTerms$X[ hits])
		bigYYes <- append( bigYYes, idTerms$Y[ hits])

		noHits <- which( cnts == 0)
		bigTileNo <- append( bigTileNo, idTerms$Tile[ noHits])
		bigXNo <- append( bigXNo, idTerms$X[ noHits])
		bigYNo <- append( bigYNo, idTerms$Y[ noHits])
	}

	close( conIn)
	cat( "\nN_reads in file:  ", nreads)

	return( list( "TileYes"=bigTileYes, "XYes"=bigXYes, "YYes"=bigYYes, "TileNo"=bigTileNo, 
			"XNo"=bigXNo, "YNo"=bigYNo))
}


`SOLiD2fastq` <- function( csfasta.file, qual.file=sub( ".csfasta", "_QV.qual", csfasta.file), 
			 	fq.file=sub( ".csfasta", ".fq.gz", csfasta.file), comment.lines=3) {

	conIn <- file( csfasta.file, open="rt")
	cat( "\nReading colorspace fasta file: ", csfasta.file)
	conInQ <- file( qual.file, open="rt")
	cat( "\nReading colorspace quality file: ", qual.file)
	conOut <- gzfile( fq.file, open="w")

	# get past the comments...
	for(i in 1:comment.lines) {
		id1 <- readLines( conIn, n=1)
		id2 <- readLines( conInQ, n=1)
		if ( substr( id1, 1,1) != "#") break
	}

	chunkSize = 200000
	cat( "\nConverting..")
	nout <- 0
	repeat {
		txt <- readLines( conIn, n=chunkSize)
		txtQ <- readLines( conInQ, n=chunkSize)
		if ( length( txt) < 1) break

		who <- seq( 1,length(txt), by=2)
		id1 <- txt[ who]
		id2 <- txtQ[ who]
		seq <- txt[ who + 1]
		qualstr <- txtQ[ who + 1]

		# make the id look like an early 'Rosetta' style of Illumina ID
		id <- base::sub( ">", "@1_", id1)
		#seq <- convColorSpaceToACGT( seq)
		seq <- fastCS2DNA( seq)
		qualstr <- base::sapply( qualstr, solexaToPhred)

		writeLines( base::paste( id, seq, "+", qualstr, sep="\n"), con=conOut)
		if ( id1[1] != id2[1]) cat( "\rBad IDs: ", id1[1], id2[1])
		nout <- nout + length(who)
		cat( ".")
		if ( nout %% 1000000 == 0) cat( formatC(nout, format="d", big.mark=","))
	}
	close( conIn)
	close( conInQ)
	close( conOut)
	cat( "\nWrote file: ", fq.file, "\nN_Reads: ", formatC( nout, format="d", big.mark=","), "\n")
	return(NULL)
}


`convColorSpaceToACGT` <- function( txt) {

	A_JUMP <- c( "A","C","G","T","N")
	C_JUMP <- c( "C","A","T","G","N")
	G_JUMP <- c( "G","T","A","C","N")
	T_JUMP <- c( "T","G","C","A","N")

	v <- strsplit( txt, split="")
	ans <- base::sapply( v, function(x) {
		N <- length(x)
		cur <- x[1]
		x[ x == "."] <- '4'
		xin <- as.numeric( x[2:N]) + 1
		out <- vector()
		now <- ""
		for ( j in 1:(N-1)) {
			now <- if ( cur == 'A') A_JUMP[ xin[j]]
				else if ( cur == 'C') C_JUMP[ xin[j]]
				else if ( cur == 'G') G_JUMP[ xin[j]]
				else if ( cur == 'T') T_JUMP[ xin[j]]
				else 'N'
			out[j] <- now
			if ( now != "N") cur <- now
		}
		return( base::paste( out, collapse=""))
		})

	return(ans)
}


