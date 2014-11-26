# pipe.VariantCalls.R

`pipe.VariantCalls` <- function( sampleIDset, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID="Pf3D7", results.path=NULL,
				seqIDset=NULL, prob.variant=0.75,
				mpileupArgs="", vcfArgs="", comboSamplesName="Combined", verbose=TRUE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	target <- getOptionValue( optT, "targetID", notfound="HsPf")
	setCurrentTarget( target)

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		results.path <- results.path
	}
	fastaFile <- getOptionValue( optT, "genomicFastaFile")

	if ( multicore.currentCoreCount() < 2) {
		nCores <- as.integer( getOptionValue( optT, "nCores", notfound=4))
		multicore.setup( nCores)
	}

	finalName <- sampleIDset[1]
	if (length(sampleIDset) > 1) finalName <- comboSamplesName
	vcfPath <- file.path( results.path, "VariantCalls", finalName)
	if ( ! file.exists( vcfPath)) dir.create( vcfPath, recursive=TRUE)

	bamfilelist <- vector()
	for ( sampleID in sampleIDset) {
		# make sure we have the BAM file already sorted
		bamfile <- paste( sampleID, "genomic.bam", sep=".")
		bamfile <- file.path( results.path, "align", bamfile)
		sortedbamfile <- BAM.verifySorted( bamfile, index=TRUE)
		bamfilelist <- c( bamfilelist, sortedbamfile)
	}

	if ( is.null( speciesID)) {
		speciesSet <- getCurrentTargetSpecies()
	} else { 
		speciesSet <- speciesID
	}
	allIDs <- allCounts <- vector()

	
	`variantCallOneSeq` <- function( sid) {
		ans <- BAM.variantCalls( bamfilelist, seqID=sid, fastaFile=fastaFile, 
					prob.variant=prob.variant, mpileupArgs=mpileupArgs, 
					vcfArgs=vcfArgs, verbose=verbose)
		cat( "\n", sid, "\tN_Variants: ", nrow(ans),"\n")
		N <- nrow(ans)
		outfile <- paste( finalName, sid, "VCF.txt", sep=".")
		outfile <- file.path( vcfPath, outfile)
		file.delete( outfile)
		if ( nrow(ans) > 0) write.table( ans, outfile, sep="\t", quote=FALSE, row.names=FALSE)
		return( N)
	}
	

	for ( speciesID in speciesSet) {
		setCurrentSpecies( speciesID)
		seqMap <- getCurrentSeqMap()
		# order to speed up the parallel computation
		seqMap <- seqMap[ order( seqMap$LENGTH, decreasing=TRUE), ]
		seqIDs <- seqMap$SEQ_ID
		if ( ! is.null( seqIDset)) seqIDs <- intersect( seqIDs, seqIDset)
		N <- length(seqIDs)
		if ( N < 1) next
		vCounts <- rep.int( 0, N)

		cat( "\n")
		ans <- multicore.lapply( seqIDs, FUN=variantCallOneSeq)
		vCounts <- as.integer( unlist( ans))
		
		allIDs <- c( allIDs, seqIDs)
		allCounts <- c( allCounts, vCounts)

		# if all the chromosomes were done, go ahead and summarize too
		if ( length( sampleIDset) == 1 && is.null( seqIDset)) {
			pipe.VariantSummary( sampleIDset[1], speciesID, annotationFile=annotationFile,
					optionsFile=optionsFile, results.path=results.path)
		}
	}
	if ( length( allIDs) < 1) return( NULL)

	out <- data.frame( "SEQ_ID"=allIDs, "N_VARIANTS"=allCounts, stringsAsFactors=F)
	return( out)
}


pipe.VariantSummary <- function( sampleID, speciesID="Pf3D7", annotationFile="Annotation.txt", 
				optionsFile="Options.txt", results.path=NULL,
				min.depth=1, min.score=5, exonOnly=FALSE, snpOnly=FALSE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	vcfPath <- file.path( results.path, "VariantCalls", sampleID)

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()
	seqIDs <- getCurrentSeqMap()$SEQ_ID
	exonMap <- getCurrentExonMap()

	out <- data.frame()

	for ( sid in seqIDs) {
		infile <- paste( sampleID, sid, "VCF.txt", sep=".")
		infile <- file.path( vcfPath, infile)
		if ( ! file.exists(infile)) {
			cat( "\rVariant Calls file not found:  ", infile, "  Skip..\n")
			next
		}
		tbl <- read.delim( infile, as.is=T)
		if ( ! all( tbl$FORMAT == "GT:PL:DP")) stop( "Required VCF format of GT:PL:DP not present...")
		tbl <- tbl[ , -match( c("FILTER","FORMAT"), colnames(tbl))]
		colnames(tbl)[ ncol(tbl)] <- "GENOTYPE_CALL"

		if (exonOnly) {
			emap <- subset.data.frame( exonMap, SEQ_ID == sid)
			tbl <- exonVariantsOnly( tbl, exonMap=emap)
		}
		
		if (snpOnly) {
			tbl <- snpVariantsOnly( tbl, mode="single")
		}

		depth <- VCF.TotalDepth( tbl$INFO)
		score <- VCF.Score( tbl$GENOTYPE_CALL)
		tbl$INFO <- depth
		colnames(tbl)[ match( "INFO",colnames(tbl))] <- "DEPTH"
		tbl$SCORE <- score
		drops <- which( score < min.score | depth < min.depth)
		if ( length( drops) > 0) tbl <- tbl [ -drops, ]

		out <- rbind( out, tbl)
		cat( "\r", sid, nrow(tbl))
	}
	
	outfile <- paste( sampleID, prefix, "Summary.VCF.txt", sep=".")
	outfile <- file.path( vcfPath, outfile)
	write.table( out, outfile, sep="\t", quote=F, row.names=F)
	cat( "\nWrote variant call summary: ", outfile, "\nN_Calls: ", nrow(out), "\n")

	return()
}


`exonVariantsOnly` <- function( tbl, exonMap=getCurrentExonMap()) {

	if ( nrow(tbl) < 1) return(tbl)

	mySID <- unique.default( tbl$SEQ_ID) 
	if ( length( mySID) > 1) stop( "Error in 'exonVariantsOnly':  Expected exactly one chromosome of SNP data")
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
	snpPos <- findInterval( tbl$POSITION, posVec, all.inside=FALSE)

	# now call each SNP position as being in an exon or not
	isExon <- sapply( snpPos, function(x) {
		
		# outside the boundaries is instant NO
		if ( x < 1 || x >= NE2) return( FALSE)
		# in an exon is NO
		if ( typeVec[ x] == "E") return( TRUE)
		FALSE
	})

	isIntron <- ! isExon
	if ( any( isIntron)) {
		cat( "  N_NonExon: ", sum( isIntron))
		return( tbl[ isExon, ])
	} else {
		return( tbl)
	}
}


`pipe.VariantComparison` <- function( sampleIDset, groupSet=sampleIDset, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID="Pf3D7", results.path=NULL, 
				min.deltaScore=40, exonOnly=FALSE, snpOnly=FALSE) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	vcfPath <- file.path( results.path, "VariantCalls")

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# gather all the VCF calls from all these samples
	seqV <- posV <- geneV <- refV <- altV <- scoreV <- qualV <- sidV <- gidV <- depthV <- vector()
	nall <- 0

	NS <- length(sampleIDset)
	cat( "\nLoading all VCF summary data..\n")
	for ( i in 1:NS) {
		sampleID <- sampleIDset[i]
		groupID <- groupSet[i]

		vcffile <- file.path( vcfPath, sampleID, paste( sampleID, prefix, "Summary.VCF.txt", sep="."))
		if ( ! file.exists( vcffile)) {
			cat( "VCF Summary file not found: ", vcffile)
			return(NULL)
		}

		tbl <- read.delim( vcffile, as.is=T)
		NC <- nrow(tbl)
		if ( NC < 1) next

		if (snpOnly) {
			tbl <- snpVariantsOnly( tbl, mode="single")
			NC <- nrow(tbl)
			if ( NC < 1) next
		}

		if (exonOnly) {
			# clean it by SeqID
			sids <- sort( unique( tbl$SEQ_ID))
			newtbl <- data.frame()
			for ( s in sids) {
				isS <- which( tbl$SEQ_ID == s)
				newsml <- exonVariantsOnly( tbl[ isS, ])
				newtbl <- rbind( newtbl, newsml)
			}
			tbl <- newtbl
			NC <- nrow(tbl)
			if ( NC < 1) next
		}

		now <- (nall+1) : (nall+NC)
		seqV[ now] <- tbl$SEQ_ID
		posV[ now] <- tbl$POSITION
		geneV[ now] <- tbl$GENE_ID
		refV[ now] <- tbl$REF_BASE
		altV[ now] <- tbl$ALT_BASE
		qualV[ now] <- tbl$QUAL
		depthV[ now] <- tbl$DEPTH
		scoreV[ now] <- tbl$SCORE
		sidV[ now] <- sampleID
		gidV[ now] <- groupID
		nall <- max( now)
		cat( "\n", i, sampleID, nall)
	}
	cat( "\nDone.\n")

	# any tweaking of fields...
	# 1. the ALT base may be comma separated...  keep only the first
	altV <- sub( ",.+","", altV)

	# now that we have all the VCF data, we need to visit each loci and summarize any 
	# differences found between groups
	cat( "\nOranizing all loci..")
	posKey <- paste( seqV, posV, sep=":")
	keyFac <- factor( posKey)
	nout <- nlevels( keyFac)
	groupFac <- factor( groupSet)
	NG <- nlevels(groupFac)
	groupNames <- levels(groupFac)
	gBase <- gScore <- gDepth <- vector( length=NG)

	# storage for the final answers
	seqOut <- posOut <- geneOut <- refOut <- depthOut <- deltaOut <- vector( length=nout)
	grpBaseOut <- matrix( "", nrow=nout, ncol=NG)
	grpScoreOut <- matrix( 1, nrow=nout, ncol=NG)
	grpDepthOut <- matrix( 0, nrow=nout, ncol=NG)
	colnames(grpBaseOut) <- colnames(grpScoreOut) <- colnames(grpDepthOut) <- groupNames

	cat( "\nEvaluating ", nout, " VCF loci..\n")
	iout <- 0
	tapply( 1:nall, keyFac, function(x) {

		# start from the null state that all are the reference
		myBase <- rep.int( refV[x[1]], NS)
		myScore <- rep.int( 1, NS)
		myDepth <- rep.int( 0, NS)
		# now fill in the ALTs that we saw
		who <- match( sidV[x], sampleIDset)
		myBase[who] <- altV[x]
		myScore[who] <- scoreV[x]
		myDepth[who] <- depthV[x]

		# now we combine by group to get the average of each
		igrp <- 0
		tapply( 1:NS, groupFac, function(j) {
			allBases <- myBase[j]
			allScores <- myScore[j]
			igrp <<- igrp + 1
			gScore[igrp] <<- mean( allScores)
			gBase[igrp] <<- names( sort.int( table( allBases), decreasing=T))[1]
			gDepth[igrp] <<- mean( myDepth[j])
		})
		delta <- diff( range( gScore))
		if (NG == 1) delta <- gScore[1]

		# fill the results data for this loci
		iout <<- iout + 1
		seqOut[iout] <<- seqV[ x[1]]
		posOut[iout] <<- posV[ x[1]]
		geneOut[iout] <<- geneV[ x[1]]
		refOut[iout] <<- refV[ x[1]]
		deltaOut[iout] <<- delta
		grpBaseOut[ iout, ] <<- gBase
		grpScoreOut[ iout, ] <<- gScore
		grpDepthOut[ iout, ] <<- gDepth
		if ( iout %% 1000 == 0) cat( "\r", iout, seqOut[iout], shortGeneName( geneOut[iout], keep=1), 
						posOut[iout], gBase, gScore, gDepth, "      ")
	})

	# now build the final result
	out <- data.frame( "SEQ_ID"=seqOut, "POSITION"=posOut, "GENE_ID"=geneOut, "PRODUCT"=gene2Product( geneOut),
			"REF_BASE"=refOut, "DELTA_SCORE"=round(deltaOut, digits=2), stringsAsFactors=FALSE)
	for ( j in 1:NG) {
		small <- data.frame( grpBaseOut[ ,j], round( grpScoreOut[ ,j], digits=2), round( grpDepthOut[ ,j]),
				stringsAsFactors=FALSE)
		colnames(small) <- paste( groupNames[j], c( "BASE","SCORE","DEPTH"), sep="_")
		out <- cbind( out, small, stringsAsFactors=FALSE)
	}

	# when there are more than one different SNP base at the same loci, that should get a higher score 
	# than returned by 'diff(range())'.
	cat( "\nRescoring multiple alternate alleles..")
	if ( NG > 1) {
		altcolumns <- seq( 7, ncol(out), by=3)
		for ( i in 1:nrow(out)) {
			ref <- out$REF_BASE[i]
			alts <- substr( as.character(out[ i, altcolumns]), 1,1)
			uniqueAlts <- setdiff( unique( alts), ref)
			if ( length( uniqueAlts) > 1) {
				isAlt <- which( alts %in% uniqueAlts)
				altScores <- as.numeric( out[ i, altcolumns+1])
				altFac <- factor( alts[ isAlt])
				altScores <- tapply( altScores[isAlt], altFac, mean)
				out$DELTA_SCORE[i] <- sum( altScores)
			}
		}
	}

	# final order is those most different bwtween the groups
	keep <- which( out$DELTA_SCORE >= min.deltaScore)
	out <- out[ keep, ]
	ord <- order( -(out$DELTA_SCORE), out$SEQ_ID, out$POSITION)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)
	cat( "\nDone.\n")

	return(out)
}


`snpVariantsOnly` <- function( tbl, mode=c("compare", "single"), min.depth=5) {

	mode <- match.arg( mode)
	if ( nrow(tbl) < 1) return( tbl)

	ref <- tbl$REF_BASE
	if ( mode == "compare") {
		alt <- alt1 <- tbl[[ 7]]
		depth1 <- as.numeric( tbl[[9]])
		alt2 <- tbl[[ 10]]
		depth2 <- as.numeric( tbl[[12]])
	} else {
		alt <- sub( ",.+", "", tbl$ALT_BASE)
	}

	# we want to leave out any indels, and keep only simple SNPs
	isSNPref <- which( nchar( ref) == 1)
	isSNPalt <- which( nchar( alt) == 1)
	if ( mode == "compare") {
		isSNPalt2 <- which( nchar( alt2) == 1)
		isSNPalt <- intersect( isSNPalt, isSNPalt2)
		# secondly, are these 2 alts really different from each other
		base1 <- alt1[ isSNPalt]
		base2 <- alt2[ isSNPalt]
		good <- which( base1 != base2)

		# secondly, if a sample has too few reads, it not proof of a SNP
		#deep1 <- which( depth1[isSNPalt] >= min.depth)
		#deep2 <- which( depth2[isSNPalt] >= min.depth)
		#deepEnough <- intersect( deep1, deep2)
		#good <- intersect( good, deepEnough)

		# recall that a location matching the reference was not called a SNP, so we have no depth info...!!
		mostDeep <- pmax( depth1[isSNPalt], depth2[isSNPalt])
		deepEnough <- which( mostDeep >= min.depth)
		good <- intersect( good, deepEnough)

		# OK, its a real SNP
		isSNPalt <- isSNPalt[ good]
	} else {
		# secondly, is this alt really different from reference
		base1 <- tbl$REF_BASE[ isSNPalt]
		base2 <- tbl$ALT_BASE[ isSNPalt]
		good <- which( base1 != base2)
		isSNPalt <- isSNPalt[ good]
	}

	# those that satisfy all are the ones we want
	keep <- intersect( isSNPref, isSNPalt)
	out <- tbl[ keep, ]
	if ( length(out) > 0) rownames(out) <- 1:nrow(out)
	return( out)
}


`pipe.VariantCompare2html` <- function( tbl, outfile="VariantCompare.html", 
					sampleIDset, groupSet=sampleIDset, 
					annotationFile="Annotation.txt", optionsFile="Options.txt", 
					speciesID=getCurrentSpecies(), results.path=NULL, Ngenes=200,
					tailWidth=26, indelCharLen=12, out.path=sub( ".html", "", outfile), 
					showDepth=TRUE, label="", ...) {

	REF_column <- 5
	SNP_columns <- seq( 7, ncol(tbl), by=3)
	# some SNPs will be the same call, just different depths... ignore those
	ndiff <- apply( as.matrix( tbl[ ,SNP_columns]), MARGIN=1, function(x) length( unique.default( x)))
	same <- which( ndiff == 1)
	if ( length( same) > 0) tbl <- tbl[ -same, ]
	N <- min( Ngenes, nrow( tbl))
	if ( N < 1) return()

	optT <- readOptionsTable( optionsFile)
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		results.path <- results.path
	}
	outPath <- file.path( results.path, "VariantCalls", out.path)
	if ( ! file.exists(outPath)) dir.create( outPath, recursive=T)
	localPlotPath <- paste( out.path, "pngPlots", sep=".")
	plotPath <- file.path( outPath, localPlotPath)
	if ( ! file.exists(plotPath)) dir.create( plotPath, recursive=T)

	out <- tbl[ 1:N, ]
	toplot <- out[ order( out$SEQ_ID, out$POSITION), ]

	# format up the table
	# convert the GeneIDs to have the genomic coord...
	newID <- paste( shortGeneName( out$GENE_ID, keep=1), out$POSITION, sep=".")
	out$GENE_ID <- newID
	# trim long indels
	for ( j in c( REF_column, SNP_columns)) {
		txt <- out[[j]]
		islong <- which( nchar(txt) > indelCharLen)
		if ( length(islong) > 0) {
			newtxt <- paste( substr( txt[islong], 1, indelCharLen), "...", sep="")
			out[[j]][ islong] <- newtxt
		}
	}
	# modify the SNP column names
	colnames(out)[5:ncol(out)] <- gsub( "_", " ", colnames(out)[5:ncol(out)])

	if ( ! showDepth) {
		drops <- grep( "DEPTH", colnames(out))
		if ( length(drops)) out <- out[ , -drops]
	}

	globalOutfile <- file.path( outPath, outfile)
	table2html( out, globalOutfile, title=paste( "Variant Comparison: &nbsp; ", label), maxRows=N,
			linkPaths=localPlotPath)
	# write a text version too
	globalOutfile <- file.path( outPath, sub( "html$", "txt", outfile))
	write.table( tbl, globalOutfile, sep="\t", quote=F, row.names=F)

	for ( i in 1:N) {
		myPos <- toplot$POSITION[i]
		mySeq <- toplot$SEQ_ID[i]
		myGene <- shortGeneName( toplot$GENE_ID[i], keep=1)
		pipe.PlotSNP( sampleIDset, mySeq, myPos, tailWidth=tailWidth, groups=groupSet,
				annotationFile=annotationFile, optionsFile=optionsFile, results.path=results.path,
				asPNG=TRUE, pngPath=plotPath, label=label, ...)
		cat( "\r", i, mySeq, myGene, myPos)
	}
}


VCF.Score <- function( PLterms) {

	# we 'are back to' doing just the homozygous reference value

	# now the genotype field is GT:PL:DP, so we need the middle one
	#scoreStr <- sub( ":.+", "", PLterms)
	myPLterms <- strsplit( PLterms, split=":", fixed=T)
	scoreStr <- sapply( myPLterms, function(x) x[2])
	scoreTerms <- strsplit( scoreStr, split=",", fixed=T)
	scores <- sapply( scoreTerms, function(x) {
		max( as.integer(x[1]), na.rm=T)
	})
	scores[ is.na(scores)] <- 0
	return( scores)
}


VCF.Pvalue <- function( PLterms) {

	scores <- VCF.Score( PLterms)
	pval <- 10 ^ (scores/-10)
	return( pval)
}


VCF.HighScoreDepth <- function( PLterms) {

	highScoringDepth <- sub( ".+:", "", PLterms)
	highScoringDepth <- as.integer( highScoringDepth)
	highScoringDepth[ is.na(highScoringDepth)] <- 0
	return( highScoringDepth)
}


VCF.TotalDepth <- function( info) {

	totalDepth <- sub( "(^.*)(DP=)([0-9]+)(;.*$)", "\\3", info)
	totalDepth <- as.integer( totalDepth)
	totalDepth[ is.na(totalDepth)] <- 0
	return( totalDepth)
}


`pipe.VariantMatrix` <- function( sampleIDset, annotationFile="Annotation.txt",
				optionsFile="Options.txt", speciesID="Pf3D7", results.path=NULL,
				seqIDset=NULL, min.depth=3, mpileupArgs="") {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	target <- getOptionValue( optT, "targetID", notfound="HsPf")
	setCurrentTarget( target)

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		results.path <- results.path
	}
	fastaFile <- getOptionValue( optT, "genomicFastaFile")

	if ( multicore.currentCoreCount() < 2) {
		nCores <- as.integer( getOptionValue( optT, "nCores", notfound=4))
		multicore.setup( nCores)
	}

	vcfPath <- file.path( results.path, "VariantCalls")

	bamfilelist <- vector()
	for ( sampleID in sampleIDset) {
		# make sure we have the BAM file already sorted
		bamfile <- paste( sampleID, "genomic.bam", sep=".")
		bamfile <- file.path( results.path, "align", bamfile)
		sortedbamfile <- BAM.verifySorted( bamfile, index=TRUE)
		bamfilelist <- c( bamfilelist, sortedbamfile)
	}

	allIDs <- allCounts <- vector()

	
	`variantMatrixOneSeq` <- function( sid) {
		ans <- BAM.mpileup( bamfilelist, seqID=sid, fastaFile=fastaFile, 
					min.depth=min.depth, mpileupArgs=mpileupArgs, 
					summarize.calls=FALSE)
		N <- nrow(ans)
		cat( "\n", sid, "\tN_Bases: ", N, "\n")

		# turn this table of pileup strings into one matrix of 'percent each nucleotide'...
		return(  as.variant.matrix( ans, min.depth=min.depth))
	}
	

	setCurrentSpecies( speciesID)
	seqMap <- getCurrentSeqMap()
	# order to speed up the parallel computation
	seqMap <- seqMap[ order( seqMap$LENGTH, decreasing=TRUE), ]
	seqIDs <- seqMap$SEQ_ID
	if ( ! is.null( seqIDset)) seqIDs <- intersect( seqIDs, seqIDset)
	N <- length(seqIDs)
	vCounts <- rep.int( 0, N)

	cat( "\n")
	ans <- multicore.lapply( seqIDs, FUN=variantCallOneSeq)
	vCounts <- as.integer( unlist( ans))
		
	allIDs <- c( allIDs, seqIDs)
	allCounts <- c( allCounts, vCounts)
	if ( length( allIDs) < 1) return( NULL)

	out <- data.frame( "SEQ_ID"=allIDs, "N_VARIANTS"=allCounts, stringsAsFactors=F)
	return( out)
}


