# pipe.ChIPalignment.R

# do the A to Z alignment pipeline on a sample

`pipe.ChIPalignment` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		curHost <- system( command="hostname", intern=TRUE)
		cat( "\nStarting  'ChIP Alignment Pipeline' on Sample:     ", sampleID,
			"\nSample Data Type:  \tChIP-seq",
			"\n\nHost Machine:      \t", curHost, 
			"\nStart Date/Time:    \t", date(), "\n", sep="")
	}

	target <- getOptionValue( optionsFile, "targetID", notfound="HsPf")
	setCurrentTarget( target)
	cat( "\nTarget Species Set:    ", target, "\t", getCurrentTargetSpecies(), "\n")

	resultsPath <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)

	# file to process comes from annotation and options...
	annT <- readAnnotationTable( annotationFile)
	inputFastqFile <- getAnnotationValue( annT, sampleID, "Filename")
	fastqPath <- getOptionValue( optionsFile, "fastqData.path")
	dir <- system( paste( " ls  ", fastqPath), intern=TRUE)
	cat( "\n\n'Fastq' Folder: ", fastqPath, "\nN_Files Found: ", length(dir),  "\n")
	print( head( dir))

	# allow paired end reads, but do them in one pass
	pairedEnd <- getAnnotationTrue( annT, sampleID, "PairedEnd", notfound=FALSE)
	inputFastqFileSet <- strsplit( inputFastqFile, split=", *")[[1]]
	#if (pairedEnd) {
	#	inputFastqFileSet <- strsplit( inputFastqFile, split=", *")[[1]]
	#	nFastq <- 1
	#	sampleIDset <- getSamplePairIDs( sampleID, annotationFile)
	#} else {
	#	inputFastqFileSet <- strsplit( inputFastqFile, split=", *")[[1]]
	#	nFastq <- length( inputFastqFileSet)
	#	sampleIDset <- sampleID
	#}
	inputFastqFileSet <- file.path( fastqPath, inputFastqFileSet)
	for ( i in 1:length(inputFastqFileSet)) inputFastqFileSet[i] <- allowCompressedFileName( inputFastqFileSet[i])

	gc()
	startT <- proc.time()

	nReadsIn <- NULL
	nNoHit <- 0

	# Part 1:  ribo clearing
	# either bypasss or do it, and gather the counts and filenames as we go
	doRiboClearing <- FALSE
	if ( ! doRiboClearing) {
		nRiboU <- nRiboM <- nRiboMoreK <- 0
		notRiboFile <- inputFastqFileSet
	} else {
		status1 <- pipe.RiboClear( inputFastqFileSet, sampleID, optionsFile=optionsFile, verbose=verbose, 
					rawReadCount=nReadsIn)
		nReadsIn <- status1$RawReads
		nNotRibo <- status1$NoHitReads
		nRiboU <- status1$UniqueReads
		nRiboM <-  status1$MultiReads
		notRiboFile <- paste( sampleID, "not.ribo.fastq.gz", sep=".")
	}
	nRibo <- ( nRiboU + nRiboM)

	# Part 2:   genomic -- always do
	nGenomicU <- nGenomicM <- 0
	genomicInFile <- notRiboFile
	notGenomicFile <- paste( sampleID, "not.genomic.fastq.gz", sep=".")
	status2 <- pipe.GenomicAlign( inputFastqFile=genomicInFile, sampleID, optionsFile=optionsFile, 
				asMatePairs=FALSE, verbose=verbose, rawReadCount=nReadsIn)

	# now we really know how many reads we started with...
	nReadsInGenomic <- status2$RawReads
	if ( is.null( nReadsIn)) {
		nReadsIn <- nReadsInGenomic
	}
	nGenomicU <- status2$UniqueReads
	nGenomicM <-  status2$MultiReads
	nGenomic <- (nGenomicU + nGenomicM)

	nNotRibo <- nReadsIn - nRibo
	nNotGenomic <- nNotRibo - nGenomic

	# Part 3:  splices
	# either bypasss or do it, and gather the counts and filenames as we go
	nNotSplice <- nNotGenomic
	doSplices <- FALSE
	if ( ! doSplices) {
		nSpliceU <- nSpliceM <- 0
		notSpliceFile <- notGenomicFile
	} else {
		spliceInFile <- notGenomicFile
		notSpliceFile <- paste( sampleID,"not.splice.fastq.gz", sep=".")
		status3 <- pipe.SpliceAlign( inputFastqFile=spliceInFile, sampleID, 
				annotationFile=annotationFile, optionsFile=optionsFile,
				verbose=verbose, rawReadCount=nReadsIn)
		nSpliceU <- status3$UniqueReads
		nSpliceM <-  status3$MultiReads
	}
	nSplice <- ( nSpliceU + nSpliceM)
	nNotSplice <- nNotGenomic - nSplice

	# final no-hits
	nNoHit <- nNotSplice
	noHitFile <- notSpliceFile
	finalNoHitFile <- paste( sampleID,"noHits.fastq.gz", sep=".")
	finalNoHitFile <- file.path( resultsPath, "fastq", finalNoHitFile)
	if ( file.exists( noHitFile)) {
		nlines <- getFileLineCount( noHitFile, sampleID, what="lineCount")
		file.rename( noHitFile, finalNoHitFile)
		quickFileLineCountRecord( finalNoHitFile, sampleID, lineCount=nlines, readCount=(nlines/4))
	}

	# convert any BAM files that need it
	pipe.ConvertAllBAMs( sampleID, annotationFile=annotationFile, optionsFile=optionsFile, 
				rawReadCount=nReadsIn, dataType="ChIP-seq", verbose=verbose)

	# now do the cleanup...
	pipe.FileCleanup( sampleID, optionsFile=optionsFile, verbose=verbose)

	myTime <- elapsedProcTime( startT, proc.time(), N=nReadsIn)

	# local function to make the results summary
	`makeAlignSummary` <- function() {
	out <- vector()
	out <- base::append( out, base::paste( "\n\nInput File:           \t", inputFastqFile))
	out <- base::append( out, base::paste( "\nN_Raw Reads:          \t", 
			prettyNum( nReadsIn, width=12, format="d", big.mark=","),"\n"))

	out <- base::append( out, base::paste( "\nNoHits File:          \t", finalNoHitFile))
	out <- base::append( out, base::paste( "\nN_NoHit Reads:        \t", 
			prettyNum( nNoHit, width=12, format="d", big.mark=","), "\t", 
			as.percent( nNoHit, big.value=nReadsIn),"\n"))

	out <- base::append( out, base::paste( "\nN_Unique Genomic:     \t", 
			prettyNum( nGenomicU, width=12, format="d", big.mark=","), "\t", 
			as.percent( nGenomicU, big.value=nReadsIn)))
	out <- base::append( out, base::paste( "\nN_Multi Genomic:      \t", 
			prettyNum( nGenomicM, width=12, format="d", big.mark=","), "\t", 
			as.percent( nGenomicM, big.value=nReadsIn)))
	out <- base::append( out, base::paste( "\nAll Genomic Reads:    \t", 
			prettyNum( nGenomic, width=12, format="d", big.mark=","), "\t", 
			as.percent( nGenomic, big.value=nReadsIn),"\n"))

	return(out)
	}  # end of local results summary


	cat( verboseOutputDivider)
	cat( "\n\nFinished 'ChIP Alignment Pipeline' on Sample:     ", sampleID, "\n")
	cat( "\nPipeline:             \t Results:")

	alignSummary <- makeAlignSummary()
	cat( alignSummary)
	summaryFile <- file.path( resultsPath, "summary", paste( sampleID, "pipeline.Summary.txt", sep="."))
	writeLines( alignSummary, con=summaryFile, sep="")

	cat( "\n\nTiming Stats: \n")
	print( myTime)
	gc()

	out <- list( "nReadsIn"=nReadsIn, "nNoHit"=nNoHit, "nRibo"=nRibo, "nGenomic"=nGenomic, "nSplice"=nSplice)

	return( out)
}

