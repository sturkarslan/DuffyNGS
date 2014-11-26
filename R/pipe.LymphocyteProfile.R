# pipe.LymphocyteProfile.R

# turn RNA-seq data into a profile of lymphocyte expression

`pipe.LymphocyteProfile` <- function( sampleID=NULL, annotationFile="Annotation.txt", optionsFile="Options.txt",
				results.path=NULL, kmerSize=73, doVelvet=TRUE, doIGBlast=TRUE, verbose=TRUE) {

	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)

	# regions for IG genes
	IGseqids <- c( "Hs_grc_02", "Hs_grc_22", "Hs_grc_14")
	IGstarts <- c( 89150000, 22370000, 106030000)
	IGstops <-  c( 90280000, 23270000, 107290000)

	# regions for TR genes
	# TRD@ is fully inside TRA@, so no explicitly given
	TRseqids <- c( "Hs_grc_14", "Hs_grc_07", "Hs_grc_07")
	TRstarts <- c( 22080000, 141990000, 38270000)
	TRstops <-  c( 23030000, 142520000, 38410000)

	# path for all results
	velvet.path <- file.path( results.path, "VelvetContigs", sampleID, "Lymphocytes")
	contigFile <- file.path( velvet.path, "contigs.fa")

	didVelvet <- FALSE
	if ( doVelvet || ! file.exists(contigFile)) {

		# step 1:   gather all aligned reads that land in/near the lymphocytes loci
		velvetFiles <- "noHits"
		cat( "\nGathering B cell aligned reads..\n")
		pipe.GatherRegionAlignments( sampleID, IGseqids, IGstarts, IGstops, asFASTQ=TRUE, fastq.keyword="IgHits")
		velvetFiles <- c( "IgHits", velvetFiles)
		cat( "\nGathering T cell aligned reads..\n")
		pipe.GatherRegionAlignments( sampleID, TRseqids, TRstarts, TRstops, asFASTQ=TRUE, fastq.keyword="TrHits")
		velvetFiles <- c( "TrHits", velvetFiles)
	
		# step 2:  create contigs of anything that may be lymphocyte reads
		pipe.VelvetContigs( sampleID, kmerSize=kmerSize, velvetFiles, folderName="Lymphocytes", makePep=T)
		didVelvet <- TRUE
	}

	# step 3:  Throw those contigs at IgBlast
	blastOutFile <- file.path( velvet.path, "IgBlast.IG.Out.txt")
	if (didVelvet || doIGBlast) {
		cat( "\n\nSearching Contigs for B cell constructs..")
		callIgBlast( fastafile=contigFile, outfile=blastOutFile, db="IG", organism="human", path="~/IgBlast", outfmt=7)
	}
	tbl <- readIgBlastOutput( infile=blastOutFile, m=7, min.blast.score=100)
	summaryFile <- file.path( velvet.path, "IgBlast.IG.Summary.txt")
	write.table( tbl, summaryFile, sep="\t", quote=F, row.names=F)
	cat( "\nWrote IG Summary file: ", summaryFile)
	ans <- profileIgBlastCoverage( tbl, min.blast.score=100)
	profile <- ans$profile
	pctIdent <- ans$identity
	colnames(pctIdent) <- paste( "Ident", colnames(pctIdent), sep="_")
	out <- cbind( round(profile, digits=4), round(pctIdent, digits=4))
	profileFile <- file.path( velvet.path, "IgBlast.IG.Profile.txt")
	write.table( out, profileFile, sep="\t", quote=F, row.names=T)
	cat( "\nWrote IG Profile file: ", profileFile)

	blastOutFile <- file.path( velvet.path, "IgBlast.TCR.Out.txt")
	if (didVelvet || doIGBlast) {
		cat( "\n\nSearching Contigs for T cell constructs..")
		callIgBlast( fastafile=contigFile, outfile=blastOutFile, db="TR", organism="human", path="~/IgBlast", outfmt=7)
	}
	tbl <- readIgBlastOutput( infile=blastOutFile, m=7, min.blast.score=100)
	summaryFile <- file.path( velvet.path, "IgBlast.TCR.Summary.txt")
	write.table( tbl, summaryFile, sep="\t", quote=F, row.names=F)
	cat( "\nWrote TCR Summary file: ", summaryFile)
	ans <- profileIgBlastCoverage( tbl, min.blast.score=100)
	profile <- ans$profile
	pctIdent <- ans$identity
	colnames(pctIdent) <- paste( "Ident", colnames(pctIdent), sep="_")
	out <- cbind( round(profile, digits=4), round(pctIdent, digits=4))
	profileFile <- file.path( velvet.path, "IgBlast.TCR.Profile.txt")
	write.table( out, profileFile, sep="\t", quote=F, row.names=T)
	cat( "\nWrote TCR Profile file: ", profileFile)

}


`pipe.PlotLymphocyteProfile` <- function( sampleID, type=c("IG", "TCR"), min.read.pct=0.5, label="",
					annotationFile="Annotation.txt", optionsFile="Options.txt",
					results.path=NULL) {

	type <- match.arg( type)
	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	velvet.path <- file.path( results.path, "VelvetContigs", sampleID, "Lymphocytes")
	profileFile <- file.path( velvet.path, paste( "IgBlast", type, "Profile.txt", sep="."))

	# the data is stored as a cbind of the profile and then the germline percentages 
	tbl <- read.delim( profileFile, as.is=T)
	N2 <- ncol(tbl)
	N <- as.integer( N2/2)
	profile <- as.matrix( tbl[ ,1:N, drop=FALSE])
	pctIdent <- as.matrix( tbl[ ,(N+1):N2, drop=FALSE])
	ans <- list( "profile"=profile, "identity"=pctIdent)

	plotIgBlastProfile( ans, sampleID=sampleID, label=label, min.read.pct=min.read.pct)
}


`pipe.PlotLymphocyteDiffProfile` <- function( sampleID1, sampleID2, type=c("IG", "TCR"), min.read.pct=0.5, label="",
					annotationFile="Annotation.txt", optionsFile="Options.txt",
					results.path=NULL) {

	type <- match.arg( type)
	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)
	velvet.path <- file.path( results.path, "VelvetContigs", sampleID1, "Lymphocytes")
	profileFile1 <- file.path( velvet.path, paste( "IgBlast", type, "Profile.txt", sep="."))
	velvet.path <- file.path( results.path, "VelvetContigs", sampleID2, "Lymphocytes")
	profileFile2 <- file.path( velvet.path, paste( "IgBlast", type, "Profile.txt", sep="."))

	# the data is stored as a cbind of the profile and then the germline percentages 
	tbl <- read.delim( profileFile1, as.is=T)
	N2 <- ncol(tbl)
	N <- as.integer( N2/2)
	profile <- as.matrix( tbl[ ,1:N, drop=FALSE])
	pctIdent <- as.matrix( tbl[ ,(N+1):N2, drop=FALSE])
	ans1 <- list( "profile"=profile, "identity"=pctIdent)
	tbl <- read.delim( profileFile2, as.is=T)
	N2 <- ncol(tbl)
	N <- as.integer( N2/2)
	profile <- as.matrix( tbl[ ,1:N, drop=FALSE])
	pctIdent <- as.matrix( tbl[ ,(N+1):N2, drop=FALSE])
	ans2 <- list( "profile"=profile, "identity"=pctIdent)

	ans <- diffIgBlastProfiles( ans1, ans2)
	plotIgBlastDiffProfile( ans, sampleID1=sampleID1, sampleID2=sampleID2, label=label, min.read.pct=min.read.pct)
}
