# pipe.DNAalignment.R

# do the A to Z alignment pipeline on a DNA-seq sample

`pipe.DNAalignment` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				convertBAMs=TRUE, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		curHost <- system( command="hostname", intern=TRUE)
		cat( "\nStarting  'DNA Alignment Pipeline' on Sample:     ", sampleID,
			"\nSample Data Type:  \tDNA-seq",
			"\n\nHost Machine:      \t", curHost, 
			"\nStart Date/Time:   \t", date(), "\n")
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

	# allow paired end reads
	pairedEnd <- getAnnotationTrue( annT, sampleID, "PairedEnd", notfound=FALSE)
	if (pairedEnd) {
		inputFastqFileSet <- strsplit( inputFastqFile, split=", *")[[1]]
		nFastq <- 2
	} else {
		inputFastqFileSet <- strsplit( inputFastqFile, split=", *")[[1]]
		nFastq <- length( inputFastqFileSet)
	}
	inputFastqFileSet <- file.path( fastqPath, inputFastqFileSet)
	for ( i in 1:nFastq) {
		inputFastqFileSet[i] <- allowCompressedFileName( inputFastqFileSet[i])
	}

	gc()
	startT <- proc.time()

	nReadsIn <- NULL
	nNoHit <- 0

	# Part 1:  ribo clearing
	doRiboClearing <- FALSE
	if ( ! doRiboClearing) {
		nRiboU <- nRiboM <- 0
		notRiboFileSet <- inputFastqFileSet
	}
	nRibo <- ( nRiboU + nRiboM)

	# Part 2:   genomic -- always do
	nGenomicU <- nGenomicM <- nGenomicMoreK <- 0
	genomicInFileSet <- notRiboFileSet
	notGenomicFile <- paste( sampleID, "not.genomic.fastq.gz", sep=".")
	status <- pipe.GenomicAlign( inputFastqFile=genomicInFileSet, sampleID, optionsFile=optionsFile, 
				asMatePairs=(nFastq==2), verbose=verbose, rawReadCount=nReadsIn)

	# now we really know how many reads we started with...
	nReadsInGenomic <- status$RawReads
	if ( is.null( nReadsIn)) {
		nReadsIn <- nReadsInGenomic
	}
	nGenomicU <- status$UniqueReads
	nGenomicM <-  status$MultiReads
	nGenomic <- (nGenomicU + nGenomicM)

	nNotRibo <- nReadsIn - nRibo
	nNotGenomic <- nNotRibo - nGenomic

	# Part 3:  splices
	nNotSplice <- nNotGenomic
	doSplices <- FALSE
	if ( ! doSplices) {
		nSpliceU <- nSpliceM <- 0
		notSpliceFile <- notGenomicFile
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
	} else if (pairedEnd) {
		# it may be paired end data that made paired end nohits
		noHitFilePair <- paste( sub( ".gz", "", noHitFile), 1:2, "gz", sep=".")
		if ( all( file.exists( noHitFilePair))) {
			# 2 .gz files in the current folder, combine and re-zip to final location
		
			cmdline <- paste( "zcat ", paste( noHitFilePair, collapse="  "), " | gzip > ", finalNoHitFile)
			cat( "\nCombining and re-zipping paired end 'NoHits' files..")
			system( cmdline)
			cat( "  Done.\n")
			nlines1 <- getFileLineCount( noHitFilePair[1], sampleID, what="lineCount")
			nlines2 <- getFileLineCount( noHitFilePair[2], sampleID, what="lineCount")
			quickFileLineCountRecord( finalNoHitFile, sampleID, lineCount=(nlines1+nlines2), readCount=((nlines1+nlines2)/4))
			file.delete( noHitFilePair)
		}
	}

	# convert the BAM files that need it
	if (convertBAMs) pipe.ConvertAllBAMs( sampleID, annotationFile=annotationFile, optionsFile=optionsFile, 
				rawReadCount=nReadsIn, dataType="DNA-seq", verbose=verbose)

	# now do the cleanup...
	pipe.FileCleanup( sampleID, optionsFile=optionsFile, verbose=verbose)

	myTime <- elapsedProcTime( startT, proc.time(), N=nReadsIn)


	# local function to make the results summary
	`makeAlignSummary` <- function() {
	out <- vector()
	out <- base::append( out, "\n")
	out <- base::append( out, base::paste( "\nInput File:           \t", inputFastqFile))
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
	cat( "\n\nFinished 'Alignment Pipeline' on Sample:     ", sampleID, "\n")
	cat( "\nPipeline:             \t Results:")

	alignSummary <- makeAlignSummary()
	cat( alignSummary)
	summaryFile <- file.path( resultsPath, "summary", paste( sampleID, "pipeline.Summary.txt", sep="."))
	writeLines( alignSummary, con=summaryFile, sep="")

	cat( "\n\nTiming Stats: \n")
	print( myTime)
	gc()

	out <- list( "nReadsIn"=nReadsIn, "nNoHit"=nNoHit, "nGenomic"=nGenomic)

	return( out)
}

