# pipe.PlotGene.R

`pipe.PlotGene` <- function( sampleIDs, genes, annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL, targetID=NULL,
				colorColumn="Color", PLOT.FUN=NULL, asPNG=FALSE, path=".", 
				keepShortGeneName=NULL, ...) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	if( is.null( targetID)) targetID <- getOptionValue( optT, "targetID", notfound="HsPf", verbose=F)
	setCurrentTarget( targetID)
	speciesSet <- getCurrentTargetSpecies()
	curSpecies <- getCurrentSpecies()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	if (asPNG) {
		if ( ! file.exists( path)) dir.create( path, recursive=T, showWarnings=F)
	}

	# turn the list of genes into the list of species to visit
	NG <- length( genes)
	gspecies <- rep.int( NA, NG)
	gptr <- rep.int( 0, NG)
	for ( speciesID in speciesSet) {
		setCurrentSpecies( speciesID)
		gmap <- getCurrentGeneMap()
		where <- match( genes, gmap$GENE_ID, nomatch=0)
		gspecies[ where > 0] <- speciesID
		gptr[ where > 0] <- where[ where > 0]

		where <- match( genes, gmap$NAME, nomatch=0)
		gspecies[ where > 0] <- speciesID
		gptr[ where > 0] <- where[ where > 0]
	}

	# at this point, we should know most all of them
	notfound <- which( is.na( gspecies))
	if ( length(notfound) > 0) {
		cat( "\nSome Gene names not found: ", genes[ notfound])
	}

	# make sure we know all these sampleIDs
	annT <- readAnnotationTable( annotationFile)
	doSamples <- intersect( sampleIDs, annT$SampleID)
	if ( length( doSamples) < 1) {
		cat( "\nNo Samples found that match annotation: ", annT$SampleID)
		return()
	}
	if ( length( doSamples) < length( sampleIDs)) {
		cat( "\nSome Samples not found: ", setdiff( sampleIDs, doSamples))
	}
	doMultiWIG <- ( length( doSamples) > 1)
	wigColors <- annT[[ colorColumn]][ match( doSamples, annT$SampleID)]

	# visit those species...
	doSpecies <- sort( unique( gspecies[ ! is.na(gspecies)]))
	for (speciesID in doSpecies) {

		# get the WIG data we need
		setCurrentSpecies(speciesID)
		gmap <- getCurrentGeneMap()
		prefix <- getCurrentSpeciesFilePrefix()
		if ( length( doSpecies) > 1) cat( "\n", speciesID, ": \n", sep="")

		if (doMultiWIG) {
			WIGlist <- loadMultipleWIGs( doSamples, speciesID, results.path)
		} else {
			wigfile <- paste( doSamples, prefix, "WIG.rda", sep=".")
			wigfile <- file.path( results.path, "wig", wigfile)
			load( wigfile)
		}

		# now visit each gene in this species
		# use chromosomal order for speed
		mygenes <- which( gspecies == speciesID)
		mygenes <- mygenes[ order( gptr[mygenes])]

		for ( ig in mygenes) {
			g <- gmap$GENE_ID[ gptr[ ig]]
			gname <- g
			if ( ! is.null( keepShortGeneName)) gname <- shortGeneName( g, keep=keepShortGeneName)

			if (doMultiWIG) {
				plotMultiWIGgene( WIGlist, colors=wigColors, gene=g, PLOT.FUN=PLOT.FUN, ...)
			} else {
				plotWIGgene( wiggles, gene=g, ...)
			}

			if ( asPNG) {
				plotfile <- paste( gname, "png", sep=".")
				plotfile <- file.cleanSpecialCharactersFromFileName( plotfile)
				plotfile <- file.path( path, plotfile)
				dev.print( png, plotfile, width=1000, height=700, "bg"="white")
			}
			if (length(mygenes) > 1) cat( "  ",gname)
		}
	}
	if ( NG > 1) cat( "\n")
		
	return( )
}

