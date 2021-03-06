# pipe.VelvetContigs.R


`pipe.VelvetContigs` <- function( sampleID, kmerSize, fastqSource=NULL,
		annotationFile="Annotation.txt",
		optionsFile="Options.txt", velvetPath="~/Velvet/velvet_1.2.10",
		minLength=200, minCoverage=3, buildHash=TRUE, buildContigs=TRUE, 
		advancedContigArgs="", folderName=NULL, 
		kmer.subfolder=FALSE, makePeptides=FALSE, verbose=TRUE) {

	require( Biostrings)

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nStarting 'VelvetContigs' for Sample:     ", sampleID, "\n\n")
	}

	optT <- readOptionsTable( optionsFile)
	resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=FALSE)

	# velvet wants a folder for writing all its results
	outpath <- file.path( resultsPath, "VelvetContigs", sampleID)
	if ( ! is.null( folderName)) {
		outpath <- file.path( resultsPath, "VelvetContigs", sampleID, folderName)
	}
	if (kmer.subfolder) outpath <- file.path( outpath, paste( "Kmer", kmerSize, sep="_"))
	if ( ! file.exists( outpath)) dir.create( outpath, recursive=T)

	# what fastq files go into the tool?
	# deafault is the full raw fastq data before any alignments are done
	if ( is.null( fastqSource)) {
		fastqPath <- getOptionValue( optT, "fastqData.path", verbose=FALSE)
		fastqFile <- getAnnotationValue( annT, key=sampleID, columnArg="Filename", verbose=FALSE)
		fastqFile <- strsplit( fastqFile, split=", *")[[1]]
		fastqFile <- file.path( fastqPath, fastqFile)
		pairedEnd <- getAnnotationTrue( annT, key=sampleID, columnArg="PairedEnd", notfound=FALSE, verbose=FALSE)
		doPairedEnd <- (pairedEnd && length(fastqFile) == 2)
	} else { 
		# otherwise, some files from the 'fastq' results subfolder for this sample
		fastqPath <- file.path( resultsPath, "fastq")
		fastqFile <- paste( sampleID, fastqSource, "fastq.gz", sep=".")
		fastqFile <- file.path( fastqPath, fastqFile)
		doPairedEnd <- FALSE
		if (verbose) {
			cat( "\nUsing Fastq files:", fastqFile, sep="\n")
		}
	}

	ans <- makeVelvetContigs( fastqFile, outpath=outpath, velvetPath=velvetPath, 
			buildHash=buildHash, buildContigs=buildContigs,
			kmerSize=kmerSize, minLength=minLength, minCoverage=minCoverage, 
			advancedContigArgs=advancedContigArgs, doPairedEnd=doPairedEnd, verbose=verbose)

	if ( makePeptides) {
		avgPepLen <- makeVelvetPeptides( sampleID, outpath=outpath, verbose=verbose)
		ans$AvgPeptide <- round(avgPepLen)
	}

	return( ans)
}


`pipe.VelvetContigSurvey` <- function( sampleID, kmerSizes, fastqSource=NULL,
		annotationFile="Annotation.txt",
		optionsFile="Options.txt", velvetPath="~/Velvet/velvet_1.2.10",
		buildHash=TRUE, buildContigs=TRUE, makePeptides=FALSE) {

	NK <- length( kmerSizes)
	if ( NK < 1) stop( "'kmerSizes' must be a vector of kmer lengths for Velvet")

	require( Biostrings)
	optT <- readOptionsTable( optionsFile)
	resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	outpath <- file.path( resultsPath, "VelvetContigs", sampleID)


	callOneSurvey <- function(kmer) {
		cat( "\n\nTrying kmer size: ", kmer)
		ans <- pipe.VelvetContigs( sampleID, kmerSize=kmer, fastqSource=fastqSource,
				annotationFile=annotationFile, optionsFile=optionsFile,
				velvetPath=velvetPath, kmer.subfolder=T, 
				buildHash=buildHash, buildContigs=buildContigs,
				makePeptides=makePeptides, verbose=FALSE)
		return( ans)
	}

	if ( length( kmerSizes) > 1) {
		multicore.ans <- multicore.lapply( kmerSizes, FUN=callOneSurvey)
		out <- data.frame()
		for ( i in 1:NK) {
			ans <- multicore.ans[[ i]]
			out <- rbind( out, ans)
		}
	} else {
		out <- callOneSurvey( kmerSizes[1])
	}

	# see if we can auto-detect the best cutoff
	out$KmerSize <- kmerSizes
	out$Rank <- 1:NK
	if ( 'AvgPeptide' %in% colnames(out)) {
		out <- out[ , c(9,10,1:8)]
	} else {
		out <- out[ , c(8,9,1:7)]
	}

	drops <- which( is.na( out$N50) | out$N50 < 10)
	if ( length( drops) > 0) out <- out[ -drops, ]

	out <- rankVelvetContigSurvey( out)

	if ( dev.cur() < 2) x11( bg="white", width=12, height=9)
	graphVelvetContigSurvey( out, outpath, label=sampleID) 
	pngfile <- file.path( outpath, paste( sampleID, "VelvetContigSurvey.png", sep="."))
	dev.print( png, pngfile, width=800, height=600)

	return( out)
}


rankVelvetContigSurvey <- function( tbl) {

	NK <- nrow( tbl)
	out <- tbl
	hasPepLen <- ( "AvgPeptide" %in% colnames(tbl))

	rankM <- matrix( 0, nrow=NK, ncol=5)
	if (hasPepLen) rankM <- matrix( 0, nrow=NK, ncol=6)

	# try to build weights from the range of the values
	rankM[ order( out$N50, decreasing=T), 1] <- 1:NK * sqrt(sd(out$N50)/mean(out$N50))
	rankM[ order( out$Max, decreasing=T), 2] <- 1:NK * sqrt(sd(out$Max)/mean(out$Max))
	rankM[ order( out$Total, decreasing=T), 3] <- 1:NK * sqrt(sd(out$Total)/mean(out$Total))
	rankM[ order( out$Usage, decreasing=T), 4] <- 1:NK * sqrt(sd(out$Usage)/mean(out$Usage))
	rankM[ order( out$Nodes, decreasing=F), 5] <- 1:NK * sqrt(sd(out$Nodes)/mean(out$Nodes))
	if (hasPepLen) {
		rankM[ order( out$AvgPeptide, decreasing=T), 6] <- 1:NK * sqrt(sd(out$AvgPeptide)/mean(out$AvgPeptide))
	}
	avgRank <- apply( rankM, 1, mean)
	out$Rank <- avgRank
	ord <- order( avgRank)
	out <- out[ ord, ]
	rownames( out) <- 1:nrow(out)
	return( out)
}


graphVelvetContigSurvey <- function( tbl, path, lwd=4, label="", legend.cex=1.2) {

	# there may have been Kmer lengths that failed, trap and remove those
	drops <- is.na( tbl$N50)
	if ( any( drops)) {
		cat( "\n\nDropping some Kmers for failing: ", tbl$KmerSize[drops])
		tbl <- tbl[ !drops, ]
	}
	if (nrow(tbl) < 1) {
		cat( "No good Kmer sizes..  Nothing to plot")
		return()
	}

	pathSet <- file.path( path, paste( "Kmer", tbl$KmerSize, sep="_"))
	bestN50 <- max( as.numeric(tbl$N50), na.rm=T)
	graphVelvetStats( path=pathSet[1], kmerSize=tbl$KmerSize[1], n50=bestN50, col=1, lwd=lwd, replot=T, label=label)

	if ( nrow(tbl) > 1) {
		for (i in 2:nrow(tbl)) graphVelvetStats( path=pathSet[i], kmerSize=tbl$KmerSize[i], 
				n50=tbl$N50[i], col=i, lwd=lwd, replot=F)
	}
	# redraw the best so its on top
	graphVelvetStats( path=pathSet[1], kmerSize=tbl$KmerSize[1], n50=tbl$N50[1], col=1, lwd=lwd, replot=F)

	legend( "topright", paste( "Kmer =", tbl$KmerSize), lwd=4, col=1:nrow(tbl), cex=legend.cex, bg="white")
}


cleanupVelvetFiles <- function( path="results/VelvetContigs") {

	# delete any large unneeded files left behind by Velvet...
	filePatterns <- c( "^Roadmaps$", "^Sequences$", "^Graph$", "^Graph2$", "^LastGraph$", "PreGraph")

	totalBytes <- 0
	nfiles <- 0

	for (patt in filePatterns) {
		files <- dir( path, pattern=patt, recursive=T, full.name=T)
		if ( length( files) < 1) next
		for ( f in files) {
			bytes <- file.info( f)$size
			nfiles <- nfiles + 1
			cat( "\r", nfiles, "  Size: ", bytes, "  File: ", f)
			totalBytes <- totalBytes + bytes
			file.delete( f)
		}
		cat( "\n")
	}
	cat( "\nDeleted_Files: ", nfiles, "\tDeleted Bytes: ", prettyNum( totalBytes, big.mark=","), "\n")
}
