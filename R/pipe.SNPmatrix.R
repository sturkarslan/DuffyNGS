# pipe.SNPmatrix.R


`pipe.SNP.BaseDepth` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				results.path=NULL, seqID=NULL, SNPtablePath="~/SNPs/") {

	BASES <- c("A","C","G","T")
	N_BASES <- 4

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	curSpecies <- getCurrentSpecies()
	geneMap <- getCurrentGeneMap()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	fastaFile <- getOptionValue( optT, "genomicFastaFile", verbose=F)

	if ( ! is.null(SNPtablePath)) SNP_curSNPpath <<- SNPtablePath

	bamfile <- file.path( results.path, "align", paste( sampleID, "genomic.bam", sep="."))
	bamfile <- BAM.verifySorted( bamfile)

	snpfolder <- file.path( results.path, "VariantCalls", sampleID)
	if ( ! file.exists( snpfolder)) dir.create( snpfolder, recursive=T)
	snpfile <- file.path( snpfolder, paste( sampleID, "SNP.BaseDepth.txt", sep="."))

	allSeqID <- getCurrentSeqMap()$SEQ_ID
	if ( ! is.null( seqID)) allSeqID <- seqID


	getSNPmatrixOneSeq <- function( seqID) {

		# gmap <- subset.data.frame( geneMap, REAL_G == TRUE & SEQ_ID %in% allSeqID)
		gmap <- subset.data.frame( geneMap, REAL_G == TRUE & SEQ_ID == seqID)
		allGenes <- gmap$GENE_ID

		# we can pre-build the storage for the answer
		snpSeq <- snpPos <- vector()
		loadKnownSNPtable( seqID)
		if (nrow(SNP_curSNPtable) < 1) return(NULL)
		keepers <- which( SNP_curSNPtable$GENE_ID %in% allGenes)
		snpSeq <- SNP_curSNPtable$SEQ_ID[ keepers]
		snpPos <- SNP_curSNPtable$POSITION[ keepers]
		snpGene <- SNP_curSNPtable$GENE_ID[ keepers]
		nSNP <- length( snpPos)	

		cat( "\nMeasuring Base Depth @ PlasmoDB SNPs:  ", seqID, "\tN_SNP: ", nSNP, "\n")
		out <- matrix( NA, nrow=N_BASES, ncol=nSNP)
		rownames(out) <- BASES
		colnames(out) <- paste( seqID, snpGene, snpPos, sep=":")

		# visit every gene, and build it up
		for ( i in 1:nrow(gmap)) {
	
			geneID <- gmap$GENE_ID[i]
			# if no SNPs to a gene, skip it now
			if ( is.na( match( geneID, SNP_curSNPtable$GENE_ID))) next

			# grab the  reads in this gene
			xLo <- gmap$POSITION[i]
			xHi <- gmap$END[i]
			curMPU <- BAM.mpileup( bamfile, seqID, fastaFile, start=xLo, stop=xHi, summarize.calls=FALSE,
					verbose=FALSE)
			if (nrow( curMPU) < 1) next

			# only keep the part at PlasmoDB defined SNPs
			curMPU <- subset.data.frame( curMPU, POSITION %in% snpPos)
			if (nrow( curMPU) < 1) next

			# turn these to the tabular form
			calls <- MPU.callBases( curMPU$CALL_BASE, curMPU$REF_BASE)
			callStrs <- MPU.callTableToString( calls$depth.table)
			matrx <- MPU.callStringsToMatrix( callStrs)
			baseCounts <- MPU.callMatrixToBaseCounts( matrx, curMPU$REF_BASE, normalize=FALSE)

			# now we can stash this data where it goes
			theseSNPs <- match( curMPU$POSITION, snpPos)
			out[ , theseSNPs] <- t( baseCounts )
			if ( i %% 100 == 0) cat( "\r", i, geneID, substr( gmap$PRODUCT[i],1,50), 
				"\tN_SNP: ", nrow(curMPU), "   ")
		}

		# send back the bases as the column
		t( out)
	}


	# do each seqID
	ans <- multicore.lapply( rev(allSeqID), getSNPmatrixOneSeq)
	names(ans) <- rev(allSeqID)

	# repackage the results
	out <- NULL
	cat( "\nCombining by chromosome..\n")
	for ( i in rev( 1:length( allSeqID))) {
		thisM <- ans[[i]]
		if ( is.null( thisM)) next

		if ( is.null( out)) {
			out <- thisM
		} else {
			out <- rbind( out, thisM)
		}
		cat( "\n", names(ans)[i], nrow(out))
	}

	# write it out
	write.table( out, snpfile, sep="\t", quote=F, row.names=TRUE)
	cat( "\nWrote SNP file:  ", snpfile)
}


`exonSNPsOnly` <- function( tbl, seqID=tbl$SEQ_ID, pos=tbl$POSITION, exonMap=getCurrentExonMap()) {

	if ( nrow(tbl) < 1) return(tbl)

	mySID <- unique.default( seqID) 
	if ( length( mySID) > 1) stop( "Error in 'exonSNPsOnly':  Expected exactly one chromosome of SNP data")
	eSID <- unique.default( exonMap$SEQ_ID) 
	if ( length( eSID) > 1) {
		exonMap <- subset.data.frame( exonMap, SEQ_ID == mySID)
		if ( nrow( exonMap) < 2) return(tbl)
	}
	
	# build a lookup table of exon starts and stops
	NE2 <- nrow( exonMap) * 2
	posVec <- rep.int( 0, NE2)
	geneVec <- rep.int( "", NE2)
	typeVec <- rep.int( "", NE2)
	i <- seq.int( 1, NE2, by=2)
	posVec[ i] <- exonMap$POSITION
	posVec[ i+1] <- exonMap$END      #+ 1
	typeVec[ i] <- "E"
	typeVec[ i+1] <- "I"
	geneVec[ i] <- exonMap$GENE_ID
	geneVec[ i+1] <- exonMap$GENE_ID

	# force to increasing order for 'findInterval'
	ord <- order( posVec)
	posVec <- posVec[ ord]
	typeVec <- typeVec[ ord]
	geneVec <- geneVec[ ord]

	# find the location of all the SNPs
	snpPos <- findInterval( pos, posVec, all.inside=FALSE)

	# now call each SNP position as being in an exon or not
	isExon <- sapply( snpPos, function(x) {
		
		# outside the boundaries is instant NO
		if ( x < 1 || x >= NE2) return( FALSE)
		# in an exon is YES
		if ( typeVec[ x] == "E") return( TRUE)
		FALSE
	})

	isIntron <- ! isExon
	if ( any( isIntron)) {
		return( tbl[ isExon, ])
	} else {
		return( tbl)
	}
}


`pipe.SNP.FreqMatrix` <- function( sampleIDset, annotationFile="Annotation.txt", optionsFile="Options.txt", 
					results.path=NULL, na.rm=TRUE, zero.rm=TRUE, min.diff=NULL) {

	BASES <- c("A","C","G","T")
	N_BASES <- 4

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	curSpecies <- getCurrentSpecies()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	snpPath <- file.path( results.path, "VariantCalls")

	snpfiles <- file.path( snpPath, sampleIDset, paste( sampleIDset, "SNP.BaseDepth.txt", sep="."))
	finfo <- file.info( snpfiles)
	missing <- which( is.na( finfo$size))
	if ( length(missing)) {
		cat( "\nSome files not found: ", sampleIDset[missing],  "  Dropping..")
		snpfiles <- snpfiles[ -missing]
		sampleIDset <- sampleIDset[ -missing]
	}
	nFiles <- length( snpfiles)

	# get each file, and accumulate the data
	out <- NULL

	cat( "\nLoading..")
	for ( i in 1:nFiles) {

		cat( " ", sampleIDset[i])
		# read in the file of base counts
		tbl <- read.delim( snpfiles[i], as.is=T)
		if ( ! all( colnames(tbl) == BASES)) {
			cat( "\nUnexpected column names in file..  Skipping: ", basename(snpfiles[i]))
			next
		}
		tbl <- as.matrix( tbl)

		# convert to frequency
		rsum <- apply( tbl, 1, sum)
		for ( j in 1:N_BASES) tbl[ ,j] <- tbl[ ,j] * 100 / rsum

		# transform to a row for each base
		v <- as.vector( t( tbl))
		names(v) <- paste( rep( rownames(tbl), each=N_BASES), rep( BASES, times=nrow(tbl)), sep=":")
		nSNP <- length(v)

		# stash the data, either making new storage or adding to 
		if ( is.null( out)) {
			out <- matrix( NA, nrow=nSNP, ncol=nFiles)
			colnames(out) <- sampleIDset
			rownames(out) <- names(v)
		} else {
			if ( nSNP != nrow(out)) {
				cat( "\nSNP Data is wrong size!   File: ", basename(snpfiles[i]))
				next
			}
			if ( any( names(v) != rownames(out))) {
				bad <- which( names(v) != rownames(out))
				cat( "\nSNP names do not match!  ", length(bad), " of ", nrow(out), "\n")
				toShow <- data.frame( "Current"=rownames(out)[bad], "Now"=names(v)[bad])
				print( head( toShow), min( 12, length(bad)))
				next
			}
		}

		# data is OK, add it
		out[ , i] <- v
	}

	# we have it all

	# drop any empty columns
	drops <- vector()
	for ( i in 1:ncol(out)) if ( all( is.na( out[ ,i]))) drops <- c( drops, i)
	if ( length(drops)) out <- out[ , -drops, drop=FALSE]

	# drop any empty rows
	if (na.rm) {
		rsums <- apply( out, 1, sum)
		drops <- which( is.na( rsums))
		if ( length(drops)) {
			out <- out[ -drops, , drop=FALSE]
			cat( "\nDropped rows with any NA:  ", length(drops))
		}
	}

	# drop any zero rows
	if (zero.rm) {
		rmaxs <- apply( out, 1, max, na.rm=T)
		drops <- which( rmaxs == 0)
		if ( length(drops)) {
			out <- out[ -drops, , drop=FALSE]
			cat( "\nDropped rows of all Zero:  ", length(drops))
		}
	}

	# drop any minimal difference rows
	if ( ! is.null( min.diff)) {
		rdiffs <- apply( out, 1, function(x) diff( range(x)))
		drops <- which( rdiffs < min.diff)
		if ( length(drops)) {
			out <- out[ -drops, , drop=FALSE]
			cat( "\nDropped rows < min.diff:   ", length(drops))
		}
	}

	# we are realy done
	return( out)
}


`SNP.FreqMatrix2Calls` <- function( tbl) {

	BASES <- c("A","C","G","T")
	N_BASES <- 4

	snpID <- rownames(tbl)
	snpLoc <- sub( ":[ACGT]$", "", snpID)
	snpBase <- sub( ".+:[0-9]+:", "", snpID)

	snpFac <- factor( snpLoc)

	out <- matrix( NA, nrow=nlevels(snpFac), ncol=ncol(tbl))
	colnames(out) <- colnames(tbl)
	rownames(out) <- levels(snpFac)
	rowOut <- 0

	tapply( 1:nrow(tbl), snpFac, function(x) {

			rowsIn <- x
			rowOut <<- rowOut + 1
			x1 <- x[1]
			for ( j in 1:ncol(tbl)) {
				cnts <- tbl[ rowsIn, j]
				if ( all( is.na( cnts))) next
				best <- which.max( cnts)
				mybase <- snpBase[ x1 + best - 1]
				out[ rowOut, j] <<- mybase
			}
			if ( rowOut %% 100 == 0) cat( "\r", rowOut, x1, snpLoc[x1])
		})

	return( out)
}


`SNP.Calls2YesNo` <- function( charM, optionsFile="Options.txt",
				fastaPath=getOptionValue( optionsFile, "genomicFastaFile")) {

	snpID <- rownames(charM)

	snpLoc <- as.integer( sub( "(.+:)([0-9]+$)", "\\2", snpID))
	snpSeq <- sub( ":.+", "", snpID)

	seqFac <- factor( snpSeq)

	out <- matrix( 0, nrow=nrow(charM), ncol=ncol(charM))
	colnames(out) <- colnames(charM)
	rownames(out) <- rownames(charM)

	tapply( 1:nrow(charM), seqFac, function(x) {

			seqID <- snpSeq[x[1]]
			fa <- getFastaSeqFromFilePath( fastaPath, seqID)

			myLocs <- snpLoc[ x]
			refBase <- strsplit( fa, split="")[[1]][ myLocs]

			for ( j in 1:ncol(charM)) {
				myBase <- charM[ x, j]
				isRef <- (refBase == myBase)
				out[ x[ isRef], j] <<- 1
			}
		})

	return( out)
}


`pipe.SNP.MOIcall` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt", 
				results.path=NULL, min.total.depth=15, min.allele.freq=0.15, 
				ignore.vargenes=TRUE, ignore.introns=TRUE) {

	BASES <- c("A","C","G","T")
	N_BASES <- 4

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	curSpecies <- getCurrentSpecies()
	geneMap <- getCurrentGeneMap()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}

	snpfile <- file.path( results.path, "VariantCalls", sampleID, paste( sampleID, "SNP.BaseDepth.txt", sep="."))
	moifile <- file.path( results.path, "VariantCalls", sampleID, paste( sampleID, "SNP.MOIcall.txt", sep="."))
	snpTbl <- as.matrix( read.delim( snpfile, as.is=T))
	cat( "\nTotal SNP sites: ", nrow(snpTbl))

	# drop any with NA no data
	drops <- which( is.na( snpTbl[ ,1]))
	if ( length(drops)) {
		snpTbl <- snpTbl[ -drops, ]
		cat( "\nDropped for no read data: ", length(drops))
	}

	# drop any with not enough reads
	totalReads <- apply( snpTbl, MARGIN=1, sum)
	drops <- which( totalReads < min.total.depth)
	if ( length(drops)) {
		snpTbl <- snpTbl[ -drops, ]
		totalReads <- totalReads[ -drops]
		cat( "\nDropped for too few reads: ", length(drops))
	}

	snpSeq <- sub( ":.+:[0-9]+$", "", rownames(snpTbl))
	snpBase <- as.integer( sub( "(^.+:)([0-9]+$)", "\\2", rownames(snpTbl)))
	snpGene <- sub( "(^.+:)(.+)(:[0-9]+$)", "\\2", rownames(snpTbl))

	# do we ignore the var genes?
	if ( ignore.vargenes) {
		vmap <- getVargeneDomainMap()
		vgenes <- sort( unique( vmap$GENE_NAME))
		drops <- which( snpGene %in% vgenes)
		if ( length(drops)) {
			snpTbl <- snpTbl[ -drops, ]
			totalReads <- totalReads[ -drops]
			snpSeq <- snpSeq[ -drops]
			snpBase <- snpBase[ -drops]
			snpGene <- snpGene[ -drops]
			cat( "\nDropped as VarGene loci:   ", length(drops))
		}
	}

	# do we ignore the introns?
	if ( ignore.introns) {
		oldN <- nrow( snpTbl)
		emap <- getCurrentExonMap()
		newTbl <- NULL
		for (sid in sort( unique( snpSeq))) {
			who <- which( snpSeq == sid)
			smlTbl <- snpTbl[ who, ]
			pos <- snpBase[ who]
			emp <- subset.data.frame( emap, SEQ_ID == sid)
			smlNew <- exonSNPsOnly( smlTbl, seqID=sid, pos=pos, exonMap=emp)
			if ( is.null( newTbl)) {
				newTbl <- smlNew
			} else {
				newTbl <- rbind( newTbl, smlNew)
			}
		}
		snpTbl <- newTbl
		snpSeq <- sub( ":[0-9]+:.+", "", rownames(snpTbl))
		snpBase <- as.integer( sub( "(^.+:)([0-9]+)(:.+)", "\\2", rownames(snpTbl)))
		totalReads <- apply( snpTbl, MARGIN=1, sum)
		snpGene <- sub( "^.+:[0-9]+:", "", rownames(snpTbl))
		cat( "\nDropped as Intron loci:   ", oldN - nrow(snpTbl))
	}

	# OK, count how many alleles are deep enough
	nDeep <- sapply( 1:nrow(snpTbl), FUN=function(x) {
			mycut <- totalReads[x] * min.allele.freq
			return( sum( snpTbl[ x, ] >= mycut))
		})

	freqTbl <- snpTbl
	for ( i in 1:nrow(snpTbl)) {
		freqTbl[ i, ] <- snpTbl[ i, ] / totalReads[i]
	}
	freqTbl <- round( freqTbl, digits=2)
	colnames(freqTbl) <- paste( "PCT", colnames(freqTbl), sep="_")

	moiSNPs <- data.frame( "SNP_ID"=rownames(snpTbl), "N_ALLELES"=nDeep, snpTbl, freqTbl, stringsAsFactors=F)
	rownames(moiSNPs) <- 1:nrow(moiSNPs)

	cat( "\nWriting MOI summary file:  ", moifile)
	write.table( moiSNPs, moifile, sep="\t", quote=F, row.names=F)
	cat( "\nN_SNPS: ", nrow(moiSNPs))

	cat( "\n\nTable of Allele Frequencies:\n")
	print( tbl <- table( moiSNPs$N_ALLELES))
	print( as.percent( tbl, big.value=nrow(moiSNPs), digits=2))

	return( moiSNPs)
}
