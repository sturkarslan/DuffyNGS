# spliceConvert.R

# translate the alignments that came from a Splice Junction Index, back into their genomic coordinate system

`spliceConvert` <-
function( filein, fileout=sub( "bam$", "converted.bam", filein), genomefile="genomic.bam", 
		spliceMapPath=bowtie2Par("IndexPath"), spliceMapPrefix="spliceMap", 
		rawReadCount=NULL, readSense="sense", readBufferSize=1000000, sampleID="", verbose=TRUE) {

	timestart <- proc.time()

	# read one existing BAM
	if ( ! file.exists( filein)) {
		cat( "\nSplice BAM file not found: ", filein)
		return(NULL)
	}
	conSplice <- bamReader( filein)
	spliceRefData <- getRefData( conSplice)

	# we also need a bam file from the genome that we are converting back too
	if ( ! file.exists( genomefile)) {
		cat( "\nGenomic BAM file not found: ", genomefile)
		return(NULL)
	}
	conGenome <- bamReader( genomefile)
	headerGenome <- getHeader( conGenome)
	genomeRefData <- getRefData( conGenome)

	# open a new BAM file to hold the fragments after resolving the splices back to the genome
	conOut <- bamWriter( headerGenome, filename=fileout)

	ans <- calcAlignSummary( mode="setup", filename=filein, rawReadCount=rawReadCount, alignPhase="Splicing")

	# set up to do a buffer at a time
	nReads <- nReadsOut <- nUniqueReads <- nMultiReads <- nBad <- 0;
	hasMore <- TRUE

	# using local function...
	my_paste <- base::paste
	my_substr <- base::substr

	# grab a buffer
	repeat {
		if ( ! hasMore) break
		cat( "\nReadBAM..")
		chunk <- getNextChunk( conSplice, n=readBufferSize, alignedOnly=TRUE)
		nNow <- size(chunk)
		if ( nNow < 1) break
		if ( nNow < readBufferSize) hasMore <- FALSE
		nReads <- nReads + nNow
		cat( "  N_Splice: ", prettyNum( as.integer( nReads), big.mark=","))

		# the SeqID field is in spliceMap notation
		spliceIndexID <- refID2seqID( refID( chunk), refData=spliceRefData)
		# gather up the splice junction mapping data for these splices
		cat( "  junctionLookup..")
		spliceMap <- spliceJunctionLookup( spliceIndexID, spliceMapPath, spliceMapPrefix);

		# some things can be done as a chunk
		spliceStart <- position( chunk)
		spliceSeq <- alignSeq( chunk)
		spliceQual <- alignQual( chunk)
		myRIDs <- readID( chunk)
		myWts <- weight.align( chunk)
		spliceIDterms <- strsplit( spliceIndexID, split="::")
		sizeTerms <- strsplit( spliceMap$SIZES, split=",")
		startTerms <- strsplit( spliceMap$POSITIONS, split=",")
		endTerms <- strsplit( spliceMap$ENDS, split=",")
		lenTerms <- nchar( spliceSeq)
		mySIDs <- spliceMap$SEQ_ID
		myNewRefIDs <- seqID2refID( mySIDs, refData=genomeRefData)

		# OK, ready to visit each splice
		cat( "  converting..")
		newChunk <- bamRange()
		newGID <- newSPID <- wtStrs <- rep( "", times=nNow*2)
		nOutNow <- 0


		for ( i in 1:nNow) {

			# for speed, try using the .Call layer directly to avoid method lookup/load
			tmpAns <- .Call( "bam_range_get_next_align", chunk@range, FALSE, FALSE, PACKAGE="DuffyNGS")
			thisAlign <- new( "bamAlign", tmpAns)

			# invalid splices get ignored
			if ( spliceMap$GENE_ID[i] == "") next
			theseTerms <- spliceIDterms[[i]]
			if ( length( theseTerms) != 3) {
				warning( paste( "Invalid spliceID: ", spliceTable$SPLICE_INDEX_ID[i]))
				nBad <- nBad + 1
				next
			}
			thisSeq <- theseTerms[1]
			thisGene <- theseTerms[2]
			thisSplice <- theseTerms[3]

			# Create the splice in terms of chromosome coordinates
			sizes <- as.numeric( sizeTerms[[i]])
			starts <- as.numeric( startTerms[[i]])
			ends <- as.numeric( endTerms[[i]])
			if ( length( starts) != 2) {
				nBad <- nBad + 1
				next
			}

			#get the splice-centric coordinates of the read
			nbases <- lenTerms[i]
			read.start = starts[1] + spliceStart[i] - 1;
			read.end = read.start + nbases - 1;
			read.frag1end <- ends[1]
			read.frag2start <- starts[2]
			loca.start <- 1
			loca.end <- nbases
			loca.frag1end <- read.frag1end - read.start + 1
			loca.frag2start <- loca.frag1end + 1
		
			# prep all the new fields we will write
			newReadID1 <- my_paste( myRIDs[i], "::splice1", sep="")
			newReadID2 <- my_paste( myRIDs[i], "::splice2", sep="")
			newPos1 <- read.start
			newPos2 <- read.frag2start
			newAlignSeq1 <- my_substr( spliceSeq[i], loca.start, loca.frag1end);
			newAlignQual1 <- my_substr( spliceQual[i], loca.start, loca.frag1end);
			newAlignSeq2 <- my_substr( spliceSeq[i], loca.frag2start, loca.end);
			newAlignQual2 <- my_substr( spliceQual[i], loca.frag2start, loca.end);

			# write these 2 new alignments, calling C directly for speed
			tmpAns <- .Call( "bam_align_modify", thisAlign@align, as.integer(myNewRefIDs[i]), 
					as.integer(newPos1), newReadID1, newAlignSeq1, newAlignQual1, PACKAGE="DuffyNGS")
			new1 <- new( "bamAlign", tmpAns)
			.Call( "bam_range_push_back", newChunk@range, new1@align, PACKAGE="DuffyNGS")

			tmpAns <- .Call( "bam_align_modify", thisAlign@align, as.integer(myNewRefIDs[i]), 
					as.integer(newPos2), newReadID2, newAlignSeq2, newAlignQual2, PACKAGE="DuffyNGS")
			new2 <- new( "bamAlign", tmpAns)
			.Call( "bam_range_push_back", newChunk@range, new2@align, PACKAGE="DuffyNGS")
			nReadsOut <- nReadsOut + 2

			newGID[nOutNow + (1:2)] <- thisGene
			newSPID[nOutNow + (1:2)] <- thisSplice
			wtStrs[nOutNow + (1:2)] <- "1"
			if ( myWts[i] < 1) wtStrs[nOutNow + (1:2)] <- formatC( myWts[i], format="f", digits=2)
			nOutNow <- nOutNow + 2
		}

		# with all the new alignments in this new range, we can slam in some new tags
		length(newGID) <- length(newSPID) <- nOutNow
		setTag( newChunk, BAMTAG_GENEID, newGID)
		setTag( newChunk, BAMTAG_SPLICEID, newSPID)
		setTag( newChunk, BAMTAG_READWEIGHT, wtStrs)

		# now we can accumulate summary facts too
		who <- seq.int( 1, nOutNow, 2)
		mySpliceIDs <- paste( shortGeneName( newGID[who], keep=1), newSPID[who], sep="::")
		mySpeciesIDs <- getSpeciesFromSeqID( mySIDs)
		ans <- calcAlignSummary( mode="addData", chunk=chunk, geneIDs=mySpliceIDs, 
						speciesIDs=mySpeciesIDs)

		cat( "  writeBAM..")
		bamSave( conOut, newChunk)
		rm( newChunk)
		gc()
	} # end of each buffer...

	bamClose( conOut)
	bamClose( conSplice)
	bamClose( conGenome)

	# tally results
	ans <- calcAlignSummary( mode="report")
	cat( "\n", ans$textSummary, "\n")

	return( c( ans, list( "Alignments"=nReadsOut, "SplicedReads"=nReads, "NonSpliceAlignments"=nBad)))
}
