# pipe.CR_Investigate.R


`dispatch.CR_Investigate` <- function( sampleID, annotationFile="Annotation.txt",
				optionsFile="Options.txt", mode="normal",
				blastIndex=getOptionValue( optionsFile, "blastIndex", notfound="nt"),
				doCR=TRUE, doBlast=TRUE, maxReads=500000, maxTime=1000, maxCycles=10, 
				ratePerCycle=1, maxCR=4000, pause=0,
				nIterations=1000, nBest=10, results.path=NULL, verbose=TRUE) {

	commandLine <- paste( "X11( width=8, height=6, xpos=30, ypos=30, bg='white'); ",
				" pipe.CR_Investigate( sampleID=\"", 
				sampleID, "\", annotationFile=\"", annotationFile, 
				"\", optionsFile=\"", optionsFile, 
				"\", blastIndex=\"", blastIndex,
				"\", mode=\"", mode, "\", doCR=", doCR, ", doBlast=", doBlast, 
				", maxReads=", as.integer(maxReads), ", maxTime=", as.integer(maxTime),
				", maxCycles=", as.integer(maxCycles), 
				", ratePerCycle=", as.numeric(ratePerCycle), 
				", maxCR=", as.integer(maxCR), 
				", pause=", as.integer(pause), 
				", nIterations=", as.integer(nIterations), ", nBest=", as.integer(nBest),
				", results.path=", if (is.null(results.path)) "NULL" else 
				paste( "\"", results.path, "\"", sep=""),
				", verbose=", verbose,
				" )", sep="")

	spawnRcmd( commandLine, Rpackage="DuffyNGS", logFile=paste( sampleID, 
			"CR_Investigate.log.txt", sep="."))

	return()
}



`pipe.CR_Investigate` <- function( sampleID, annotationFile="Annotation.txt",
		optionsFile="Options.txt", mode=c( "normal", "QuickQC", "genomic"),
		blastIndex=getOptionValue( optionsFile, "blastIndex", notfound="nt"),
		doCR=TRUE, doBlast=TRUE, maxReads=500000, maxTime=1000, maxCycles=10, 
		ratePerCycle=1, maxCR=4000, pause=0, 
		nIterations=1000, nBest=10, results.path=NULL, verbose=TRUE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nStarting 'Consensus Reads Investigate' for Sample:     ", sampleID, "\n\n")
	}
	gc()

	startT <- proc.time()

	mode <- match.arg( mode)

	# get the file of noHit reads
	if ( is.null( results.path)) {
		resultsPath <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=verbose)
	} else {
		resultsPath <- results.path
	}
	infile <- paste( sampleID, "noHits", "fastq", sep=".") 
	infile <- file.path( resultsPath, "fastq", infile)
	if ( mode == "QuickQC") {
		infile <- file.path( resultsPath, "fastq", paste( sampleID, "QuickQC.noHits.fastq", sep="."))
	}
	if ( mode == "genomic") {
		infile <- file.path( resultsPath, "fastq", paste( sampleID, "not.genomic.fastq", sep="."))
	}
	nNoHits <- getFileLineCount( infile, sampleID) / 4

	# either reload or build the Unique Short Reads data structure
	usrFile <- file.path( resultsPath, "USR", paste( sampleID, "USR.rda", sep="."))
	if ( file.exists( usrFile)) {
		cat( "\nLoading existing USR dataset:  ", usrFile)
		load( usrFile, envir=.GlobalEnv)
	} else {
		trim5 <- as.numeric( getOptionValue( optionsFile, "trim5", notfound=0))
		trim3 <- as.numeric( getOptionValue( optionsFile, "trim3", notfound=0))
		ans <- USR_setup( infile, sampleID, resultsPath, trim5=trim5, trim3=trim3, 
				Nkeep=maxReads)
		doCR <- doBlast <- TRUE
	}

	# set up folder to hold results
	crPath <- file.path( resultsPath, "CR")
	if ( ! file.exists( crPath)) dir.create( crPath, recursive=TRUE, showWarnings=FALSE)
	crPngPath <- file.path( crPath, paste( sampleID, "pngPlots", sep="."))
	if ( ! file.exists( crPngPath)) dir.create( crPngPath, showWarnings=FALSE)

	# run the automatic CR growing tool
	crFile <- file.path( crPath, paste( sampleID, "CR.rda", sep="."))
	if ( doCR || !file.exists(crFile)) {
		bestCR <- autoRunCR( nBest=nBest, nIterations=nIterations, maxTime=maxTime, 
					maxCycles=maxCycles, ratePerCycle=ratePerCycle, 
					maxCR=maxCR, pause=pause,
					contextFile=crFile, pngPath=crPngPath,
					label=sampleID)
		CRT_best <<- bestCR
		saveCRTcontext( crFile)
	} else {
		cat( "\nLoading existing CR dataset:  ", crFile)
		load( crFile, envir=.GlobalEnv)
		bestCR <- CRT_best
	}

	# see what these maight be via Blast...
	if ( doBlast) {
	blastResults <<- NULL
	if ( ! is.null( blastIndex)) {
		blastResults <<- CRblaster( sampleID, crIDs=bestCR, optionsFile=optionsFile, 
					blastIndex=blastIndex, evalue=0.1, wordsize=11)
	}}

	# turn the Blast results into a table form text file
	if ( ! is.null( blastResults)) {

	fileout <- file.path( resultsPath, "CR", paste( sampleID, "CR.BlastSummary.txt", sep="."))
	tbl <- summarizeCRblastOutput( blastResults, fileout=fileout)

	# and into clickable HTML, with plots for each
	cat( "\nMaking HTLM table...")
	htmlFile <- sub( ".txt$", ".html", fileout)
	htmlTitle <- paste( "Consensus Reads:   'no hits' reads from sample:  ", sampleID)
	htmlPngPath <- paste( sampleID, "pngPlots", sep=".")
	# try to make the consensus sequence be line wrapped and 'cut/paste'-able
	tbl$CONSENSUS_SEQUENCE <- as.text.Fasta( as.Fasta( tbl$CR_ID, tbl$CONSENSUS_SEQUENCE),
			line.width=50)
	tbl$CONSENSUS_SEQUENCE <- sub( "\n", " <br> ", tbl$CONSENSUS_SEQUENCE)
	tbl$TEXT_DESCRIPTION_OF_BEST_HITS <- gsub( " | ", " <br> ", tbl$TEXT_DESCRIPTION_OF_BEST_HITS, 
			fixed=TRUE)
	table2html( tbl, htmlFile, title=htmlTitle, linkColumnNames="CR_ID", linkPaths=htmlPngPath)
	for( i in bestCR) {
		pngFile <- file.path( crPngPath, paste( "cr_",i,".png",sep=""))
		png( pngFile, width=1000, height=700, bg="white")
		plotOneCRT( i, label=sampleID)
		dev.off()
	}
	}
	
	# done...
	USR_cleanup()
	CR_cleanup()

	gc()
	return()
}

