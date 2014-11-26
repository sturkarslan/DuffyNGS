# pipe.Transcriptome.R

# turn the Wiggle Track pipeups into gene expression results

`pipe.Transcriptome` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				speciesID=NULL, results.path=NULL, dataType=NULL,
				altGeneMap=NULL, altGeneMapLabel=NULL, loadWIG=FALSE, verbose=TRUE,
				mode="normal", byExon=FALSE)
{

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'Transcriptome' on Sample:     ", sampleID, "\t", altGeneMapLabel, "\n")
	}
	startTotalT <- proc.time()
	grandTotalReads <- 0

	optT <- readOptionsTable( optionsFile)
	target <- getOptionValue( optT, "targetID", notfound="HsPf")
	setCurrentTarget( target)
	allSpecies <- getCurrentTargetSpecies()
	if ( ! is.null(speciesID)) allSpecies <- speciesID

	if ( is.null( results.path)) {
		resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		resultsPath <- results.path
	}
	if ( ! file.exists( resultsPath)) dir.create( resultsPath, showWarnings=FALSE)

	setCurrentTarget( target)

	# other annotation facts we need...
	origSampleID <- originalSamplePairID( sampleID, annotationFile)
	sampleID <- origSampleID
	annT <- readAnnotationTable( annotationFile)

	if ( is.null(dataType)) dataType <- getAnnotationValue( annT, origSampleID, "DataType", notfound="RNA-seq")

	if ( mode == "normal" && dataType == "ChIP-seq") {
		pipe.ChIPpeaks( origSampleID, annotationFile, optionsFile, speciesID=speciesID, 
				results.path=results.path, loadWIG=loadWIG, verbose=FALSE)
		return()
	}

	strandSpecific <- getAnnotationTrue( annT, origSampleID, "StrandSpecific", notfound=TRUE)
	useBothStrands <- ! strandSpecific
	keepIntergenics <- getAnnotationTrue( annT, origSampleID, "KeepIntergenics", notfound=TRUE)

	# during 'altGeneMap' runs, skip over species that are not covered...
	doingAltGeneMap <- FALSE
	if ( ! is.null( altGeneMap)) {
		doingAltGeneMap <- TRUE
		altSpeciesSet <- unique.default( getSpeciesFromSeqID( altGeneMap$SEQ_ID))
	}

	needLoadWIG <- loadWIG

	# do each species all the way through
	for( speciesID in allSpecies) {
	
	if ( doingAltGeneMap && (!( speciesID %in% altSpeciesSet))) next

	cat( "\n\nExtracting results for Species:  ", speciesID,"\n")
	startT <- proc.time()
	gc()
	setCurrentSpecies( speciesID)
	speciesPrefix <- getCurrentSpeciesFilePrefix()

	# find and/or load our WiggleBin data structure
	fileOutWIG <- file.path( resultsPath, "wig", paste( sampleID, speciesPrefix, "WIG.rda", sep="."))

	if ( !file.exists( fileOutWIG) || needLoadWIG) {
		
		cat( "\nLoading Wiggle Track data structure from alignments...")
		pipe.AlignToWig( sampleID, annotationFile=annotationFile, optionsFile=optionsFile,
				results.path=resultsPath, dataType=dataType, mode=mode)

		#get that data back when done... object 'wiggles'
		load( file=fileOutWIG)

		if ( mode == "normal") {
			cat( "\nCalculating Wiggle Track Metrics for Gene and Strand Specificity...\n")
			binMetricsFile <- paste( sampleID, speciesPrefix, "Metrics.txt", sep=".")
			binMetricsFile <- file.path( resultsPath, "summary", binMetricsFile)
			calcWIGmetrics( wiggles, asDataFrame=FALSE, logFile=binMetricsFile)
		}
		needLoadWIG <- FALSE

	} else {
		cat( "\nLoading pre-existing Wiggles data...")
		load( file=fileOutWIG)
	}

	# we could be given a non-standard geneMap to use...
	if ( is.null( altGeneMap)) {

		# regular...
		fileOutTrans <- paste( sampleID, speciesPrefix, "Transcript.txt", sep=".")
		pathOut <- file.path( resultsPath, "transcript")
		if ( ! file.exists( pathOut)) dir.create( pathOut, showWarnings=FALSE)
		fileOutTrans <- file.path( pathOut, fileOutTrans)
		trans <- calcWigTranscriptome( wiggles, useBothStrands=useBothStrands, 
					keepIntergenics=keepIntergenics, byExon=byExon,
					fileout=fileOutTrans)

	} else {
		# alternate...
		gmap <- altGeneMap
		if ( is.character(gmap)) {
			gmap <- read.delim( file=altGeneMap, as.is=TRUE)
		}
		if ( ! all( c("GENE_ID", "SEQ_ID") %in% colnames(gmap))) 
			stop( paste( "pipe.Transcriptome: invalid alternate geneMap",
			"does not have required GENE_ID, SEQ_ID columns"))
		if ( is.null( altGeneMapLabel) || base::nchar( altGeneMapLabel) < 1) 
				stop( "pipe.Transcriptome: missing file name text term 'altGeneMapLabel' ")

		# only do it for the right species...
		if ( getSpeciesFromSeqID( gmap$SEQ_ID[1]) == speciesID) {

			cat( "\n\nDoing alternate gene map transcriptome:  ", altGeneMapLabel, "\n")
			fileOutTrans <- paste( sampleID, speciesPrefix, altGeneMapLabel, "Transcript.txt", sep=".")

			# use a sub-folder for alternate transcripts...
			pathOut <- file.path( resultsPath, "transcript", altGeneMapLabel)
			if ( ! file.exists( pathOut)) dir.create( pathOut, showWarnings=FALSE)
			fileOutTrans <- file.path( pathOut, fileOutTrans)

			trans <- calcWigTranscriptome( wiggles, geneMap=gmap, useBothStrands=useBothStrands, 
					keepIntergenics=TRUE, byExon=FALSE, fileout=fileOutTrans)
		} else {
			cat( "\nThis species not in Alternate Gene Map... Skipping. ")
		}
	}

	reads <- WIG_getTotalRawReads( wiggles)
	grandTotalReads <- grandTotalReads + reads$Unique + reads$Multi

	myTime <- elapsedProcTime( startT, proc.time(), N=(reads$Unique + reads$Multi))
	if (verbose) {
		cat( "\n\nFinished Species:  ", speciesID)
		cat( "\n\nSpecies Timing Stats: \n")
		print( myTime)
		gc()
	}

	}  # end of all speceies...

	myTime <- elapsedProcTime( startTotalT, proc.time(), N=grandTotalReads)
	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nFinished pipe.Transcriptome:  ", sampleID, "\tSpecies set:  ", allSpecies)
		cat( "\n\nSample Timing Stats: \n")
		print( myTime)
	}

	return( sampleID)
}



`dispatch.TranscriptPlusHTML` <- function( sampleID, annotationFile="Annotation.txt",
					optionsFile="Options.txt", banner="", mode="normal",
					maxReads=NULL, pause=0, results.path=NULL) {

	commandLine <- paste( "X11( width=8, height=6, xpos=20, ypos=20, bg='white'); ",
				" pipe.TranscriptPlusHTML( sampleID=\"", 
				sampleID, "\", annotationFile=\"", annotationFile, 
				"\", optionsFile=\"", optionsFile, 
				"\", banner=\"", banner, 
				"\", mode=\"", mode, "\"", 
				", maxReads=", if (is.null(maxReads)) "NULL" else as.integer(maxReads), 
				", pause=", as.integer(pause), 
				", results.path=", if (is.null(results.path)) "NULL" else 
					paste("\"",results.path,"\"",sep=""), 
				" )", sep="")

	spawnRcmd( commandLine, Rpackage="DuffyNGS", logFile=paste( sampleID, "transcript.log.txt", sep="."))

	return()
}


`pipe.TranscriptPlusHTML` <- function( sampleID, annotationFile="Annotation.txt",
				optionsFile="Options.txt", banner="", mode="normal", loadWIG=FALSE,
				maxReads=NULL, pause=0, results.path=NULL, ...) {

	# get a bit of info for setting up
	optT <- readOptionsTable( optionsFile)
	target <- getOptionValue( optT, "targetID", notfound="HsPf")
	setCurrentTarget( target)
	allSpecies <- getCurrentTargetSpecies()

	# do the transcriptomes
	pipe.Transcriptome( sampleID, annotationFile, optionsFile, speciesID=NULL,
			results.path=results.path, loadWIG=loadWIG, mode=mode)
	
	# now make some gene plots
	nGenes <- 25
	for( s in allSpecies) {
		pipe.TranscriptToHTML( sampleID, annotationsFile, optionsFile,
				speciesID=s, results.path=results.path, 
				pause=pause, N=nGenes, label=banner, ...)
	}
	return()
}
