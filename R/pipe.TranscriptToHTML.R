# pipe.TranscriptToHTML.R

# turn a sample's  transcriptome into a HTML table with plots

`pipe.TranscriptToHTML` <- function( sampleID, speciesID=NULL, annotationFile="Annotation.txt", 
				optionsFile="Options.txt", results.path=NULL,
				N=100, minYmax=10, tailWidth=1000, pause=NULL,
				altGeneMap=NULL, altGeneMapLabel=NULL, targetID=NULL,
				plotType="box", useLog=FALSE, label="", ...){

	# open a drawing device
	if ( ! (capabilities()[ "png"] )) stop( "pipe.TranscriptToHTML:  cannot open PNG device")

	cat( "\n\nMaking HTML files for Transcript:  ", sampleID, "\tspecies: ", speciesID, "\n")

	# get the target species set...
	optT <- readOptionsTable( optionsFile)
	if ( is.null( targetID)) {
		targetID <- getOptionValue( optT, "targetID", notfound="HsPf")
	}
	setCurrentTarget( targetID)

	if ( is.null( speciesID)) {
		speciesID <- getCurrentTargetSpecies()
	}

	for( species in speciesID) {

	setCurrentSpecies( species)
	speciesPrefix <- getCurrentSpeciesFilePrefix()
	cat( "\nSpeciesID: ", species)

	if ( is.null( results.path)) {
		resultsPath <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	} else {
		resultsPath <- results.path
	}

	# get the Wiggle data
	fileWIG <- paste( sampleID, speciesPrefix, "WIG.rda", sep=".")
	fileWIG <- file.path( resultsPath, "wig", fileWIG)
	who <- load( fileWIG)
	thisWIG <- get( who, envir=environment())

	# get the transcriptome
	fileTrans <- paste( sampleID, speciesPrefix, "Transcript.txt", sep=".")
	fileTrans <- file.path( resultsPath, "transcript", fileTrans)
	if ( ! is.null( altGeneMap)) {
		fileTrans <- paste( sampleID, speciesPrefix, altGeneMapLabel, "Transcript.txt", sep=".")
		fileTrans <- file.path( resultsPath, "transcript", altGeneMapLabel, fileTrans)
	}
	thisTrans <- read.delim( fileTrans, as.is=TRUE)

	# limit the HTML table and plots to the top N genes
	N <- min( N, nrow(thisTrans))
	# further limit it if there are not that many genes with non-zero reads
	nNonZero <- sum( thisTrans$READS_M > 0)
	if ( nNonZero < N) N <- nNonZero
	theseGenes <- thisTrans$GENE_ID[ 1:N]

	# make a subfolder for the plots
	pngFolder <- linkPngFolder <- paste( sampleID, speciesPrefix, "png_Plots", sep=".")
	if ( ! is.null(altGeneMap)) {
		pngFolder <- linkPngFolder <- paste( sampleID, speciesPrefix, altGeneMapLabel, "png_Plots", sep=".")
	}
	pngFolder <- file.path( resultsPath, "html", pngFolder)
	if ( ! file.exists( pngFolder)) dir.create( pngFolder, recursive=TRUE, showWarnings=FALSE)

	# OK make the HTML file
	htmlFile <- paste( sampleID, speciesPrefix, "topGenes.html", sep=".")
	if ( ! is.null(altGeneMap)) {
		htmlFile <- paste( sampleID, speciesPrefix, altGeneMapLabel, "topGenes.html", sep=".")
	}
	htmlFile <- file.path( resultsPath, "html", htmlFile)
	title <- paste( " SampleID= ", sampleID, "<br>SpeciesID= ", species, 
			"<br>", label)
			
	table2html( thisTrans[ 1:N, ], fileout=htmlFile, title=title, linkPaths=linkPngFolder)

	# nom make those individual plots
	cat( "\nMaking gene plots of top ", N, " genes.")
	makeAllWIGgenePlots( thisWIG, geneSet=theseGenes, path=pngFolder, tailWidth=tailWidth,  
			plotType=plotType, minYmax=minYmax, useLog=useLog, label=label, 
			altGeneMap=altGeneMap, pause=pause, ...)
	
	} # end of all speciesIDs...

	return()
}

