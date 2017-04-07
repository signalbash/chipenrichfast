#' A helper function to post-process peak GRanges
#'
#' Check for overlapping input regions, sort peaks, and force peak names
#'
#' @param gr A \code{GRanges} of input peaks.
#'
#' @return A \code{GRanges} that is sorted if the \code{seqinfo} is set, and has named peaks.
postprocess_peak_grs = function(gr) {
	# Check for overlapping peaks
	reduced_gr = GenomicRanges::reduce(gr)
	if(length(gr) != length(reduced_gr)) {
		warning('Some input regions overlap. It is recommended that input regions be disjoint.')
	}

	# Sort
	gr = GenomicRanges::sort(gr)

	# Create name column, or replace what's there
	GenomicRanges::mcols(gr)$name = paste0('peak:', seq_along(gr))

	return(gr)
}

#' Read files containing peaks or genomic regions
#'
#' The following formats are fully supported via their file extensions: .bed,
#' .broadPeak, .narrowPeak, .gff3, .gff2, .gff, and .bedGraph or .bdg. BED3 through
#' BED6 files are supported under the .bed extension. Files without these extensions
#' are supported under the conditions that the first 3 columns correspond to
#' chr, start, and end and that there is either no header column, or it is
#' commented out. Files may be compressed with gzip, and so might end in .narrowPeak.gz,
#' for example. For files with extension support, the \code{rtracklayer::import()}
#' function is used to read peaks, so adherence to the mentioned file formats is
#' necessary.
#'
#' NOTE: Header rows must be commented with \code{#} to be ignored. Otherwise,
#' an error may result.
#'
#' NOTE: A warning is given if any input regions overlap. In the case of enrichment
#' testing with \code{method = 'broadenrich'}, regions should be disjoint.
#'
#' Typically, this function will not be used alone, but inside \code{chipenrich()}.
#'
#' @param file_path A path to a file with input peaks/regions. See extended
#' description above for details about file support.
#' @param genome The genome build the input regions are with respect to.
#' Usually inferred from \code{chipenrich()}. Default is NA, in which case no
#' genome is assigned to the \code{GRanges} and the peaks cannot be sorted.
#'
#' @return A \code{GRanges} with \code{mcols} matching any extra columns.
#'
#' @examples
#'
#' # Example of generic .txt file with peaks
#' file = system.file('extdata', 'test_header.txt', package = 'chipenrich')
#' peaks = read_bed(file)
#'
#' # Example of BED3
#' file = system.file('extdata', 'test_assign.bed', package = 'chipenrich')
#' peaks = read_bed(file, genome = 'hg19')
#'
#' # Example of BED3 without genome
#' file = system.file('extdata', 'test_assign.bed', package = 'chipenrich')
#' peaks = read_bed(file)
#'
#' # Example of narrowPeak without genome
#' file = system.file('extdata', 'test.narrowPeak', package = 'chipenrich')
#' peaks = read_bed(file)
#'
#' # Example of gzipped broadPeak without genome
#' file = system.file('extdata', 'test.broadPeak.gz', package = 'chipenrich')
#' peaks = read_bed(file)
#'
#' # Example of gzipped gff3 Fly peaks
#' file = system.file('extdata', 'test.gff3.gz', package = 'chipenrich')
#' peaks = read_bed(file, genome = 'dm3')
#'
#' @export
read_bed = function(file_path, genome = NA) {
	message(sprintf('Reading peaks from %s', file_path))

	# Establish format and assign extraCols and format as needed
	# See https://genome.ucsc.edu/FAQ/FAQformat for format details
	if(grepl('.narrowPeak', file_path)) {
		extraCols = c('signalValue' = 'numeric', 'pValue' = 'numeric', 'qValue' = 'numeric', 'peak' = 'integer')
		format = 'bed'
	} else if (grepl('.broadPeak', file_path)) {
		extraCols = c('signalValue' = 'numeric', 'pValue' = 'numeric', 'qValue' = 'numeric')
		format = 'bed'
	} else if (grepl('.gff3', file_path) || grepl('.gff2', file_path) || grepl('.gff', file_path)) {
		extraCols = character()
		format = 'gff'
	} else if (grepl('.bedGraph', file_path) || grepl('.bdg', file_path)) {
		extraCols = character()
		format = 'bedGraph'
	} else if (grepl('.bed', file_path)) {
		# This works for BED3 through BED6
		extraCols = character()
		format = 'bed'
	} else {
		# Catch other formats
		format = 'txt'
	}

	# Read in the file using rtracklayer::import() for anything but format == 'txt'
	# Otherwise, use read.table + load_peaks() to get the GRanges we need
	if(format %in% c('gff3','gff2','gff','bedGraph')) {
		# Extension detection takes care of the above formats
		gr = rtracklayer::import(con = file_path, genome = genome)
		# Check for overlapping peaks and warn, sort, and name the peaks
		gr = postprocess_peak_grs(gr)
	} else if (format == 'bed') {
		# narrowPeak and broadPeak aren't supported with extension detection
		# so give the fixed extraCols and use the bed format.
		gr = rtracklayer::import(con = file_path, format = format, genome = genome, extraCols = extraCols)
		# Check for overlapping peaks and warn, sort, and name the peaks
		gr = postprocess_peak_grs(gr)
	} else {
		# This catches format == 'txt' (i.e. all other extensions)
		df = read.table(file_path, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
		# Rename the first three columns
		colnames(df)[1:3] = c('chr','start','end')
		# Don't keep extra columns because we don't know what they are
		# NOTE: This includes postprocess_peak_grs()
		gr = load_peaks(dframe = df, genome = genome, keep.extra.columns = FALSE)
	}

	return(gr)
}

#' Convert a BEDX+Y data.frame and into GRanges
#'
#' Given a \code{data.frame} in BEDX+Y format, use the built-in function
#' \code{GenomicRanges::makeGRangesFromDataFrame()} to convert to \code{GRanges}.
#'
#' Typically, this function will not be used alone, but inside \code{chipenrich()}.
#'
#' @param dframe A BEDX+Y style \code{data.frame}. See \code{GenomicRanges::makeGRangesFromDataFrame}
#' for acceptable column names for appropriate conversion to \code{GRanges}.
#' @param genome The genome build used to get input regions. Usually inferred from
#' \code{chipenrich()}. Default is NA, in which case no genome is assigned to the
#' \code{GRanges} and the peaks cannot be sorted.
#' @param keep.extra.columns Keep extra columns parameter from \code{GenomicRanges::makeGRangesFromDataFrame()}.
#'
#' @return A \code{GRanges} that may or may not \code{keep.extra.columns}, and
#' that may or may not be stranded, depending on whether there is strand column
#' in the \code{dframe}.
#'
#' @examples
#'
#' # Example with just chr, start, end
#' peaks_df = data.frame(
#'	chr = c('chr1','chr2','chr3'),
#'	start = c(35,74,235),
#'	end = c(46,83,421),
#'	stringsAsFactors = FALSE)
#' peaks = load_peaks(peaks_df, genome = 'hg19')
#'
#' # Example with extra columns
#' peaks_df = data.frame(
#'	chr = c('chr1','chr2','chr3'),
#'	start = c(35,74,235),
#'	end = c(46,83,421),
#'	strand = c('+','-','+'),
#'	score = c(36, 747, 13),
#'	stringsAsFactors = FALSE)
#' peaks = load_peaks(peaks_df, genome = 'hg19', keep.extra.columns = TRUE)
#'
#' @export
load_peaks = function(dframe, genome = NA, keep.extra.columns = TRUE) {
	# Use built-in function
	if(is.na(genome)) {
		gr = GenomicRanges::makeGRangesFromDataFrame(
			df = dframe,
			keep.extra.columns = keep.extra.columns)
	} else {
		gr = GenomicRanges::makeGRangesFromDataFrame(
			df = dframe,
			seqinfo = GenomeInfoDb::Seqinfo(genome = genome),
			keep.extra.columns = keep.extra.columns)
	}

	# Check for overlapping peaks and warn, sort, and name the peaks
	gr = postprocess_peak_grs(gr)

	return(gr)
}

#' Function to read custom locus definition from file
#'
#' This function reads a tab-delimited text (with a header) file that should have
#' columns 'chr', 'start', 'end', and a column named 'gene_id' (or 'geneid') with the
#' Entrez Gene ID. If a \code{supported_genomes()} is given, then a column
#' of gene symbols named 'symbol', will be added. If an unsupported genome is
#' used there are two options: 1) Have a column named 'symbols' with the gene
#' symbols in the custom locus definition, and leave \code{genome = NA}, or 2)
#' leave \code{genome = NA}, do not provide gene symbols, and NAs will be used.
#'
#' @param file_path A file path for the custom locus definition.
#' @param genome A genome from \code{supported_genomes()}, default \code{NA}.
#'
#' @return A \code{LocusDefinition} class object with slots \code{dframe},
#' \code{granges}, \code{genome.build}, and \code{organism}.
read_ldef = function(file_path, genome = NA) {

	message("Reading user-specified gene locus definitions: ", file_path)
	df = read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

	# Remove duplicated rows, not duplicated genes.
	# Some genes will exist on multiple rows because they are split by
	# other transcripts, or small nuclear RNAs, etc.
	df = unique(df)

	# Rename 'geneid' columns 'gene_id' for backwards compatibility
	if('geneid' %in% colnames(df)) {
		colnames(df)[which(colnames(df) == 'geneid')] = 'gene_id'
	}

	if(!('gene_id' %in% colnames(df))) {
		stop("Error: Custom locus definition must have column named 'gene_id'.")
	}

	if('symbol' %in% colnames(df)) {
		# Do nothing
		message('Using given symbol column...')
	} else if (genome %in% supported_genomes()) {
		message('Using orgDb package to get gene symbols...')
		eg2symbol = genome_to_orgdb(genome)
		df$symbol = eg2symbol[match(df$gene_id, eg2symbol$gene_id), 'symbol']
	} else {
		message('Setting gene symbols to NA...')
		df$symbol = NA
	}

	if(genome %in% supported_genomes()) {
		gr = GenomicRanges::makeGRangesFromDataFrame(
			df = df,
			seqinfo = GenomeInfoDb::Seqinfo(genome = genome),
			keep.extra.columns = TRUE)
	} else {
		gr = GenomicRanges::makeGRangesFromDataFrame(
			df = df,
			keep.extra.columns = TRUE)
	}

	# Create new locus definition object.
	object = new("LocusDefinition")

	object@dframe = df
	object@granges = gr

	return(object)
}

# Creates an object that mimics the GeneSet class
# from the chipenrich.data package, as in:
#    setClass("GeneSet",representation(
#      set.gene = "environment",
#      type = "character",
#      set.name = "environment",
#      all.genes = "character",
#      organism = "character",
#      dburl = "character"
#    ),prototype(
#      set.gene = new.env(parent=emptyenv()),
#      type = "",
#      set.name = new.env(parent=emptyenv()),
#      all.genes = "",
#      organism = "",
#      dburl = ""
#    ),
#      package = "chipenrich.data"
#    )

#' Function to read custom gene sets from file
#'
#' This function reads a two-columned tab-delimited text file (with header). Column
#' names are ignored, but the first column should be geneset names or IDs and the
#' second column should be Entrez Gene IDs.
#'
#' @param file_path A file path for the custom gene set.
#'
#' @return A \code{GeneSet} class object.
read_geneset = function(file_path) {

	# Read in flat file.
	message("Reading user-specified gene set definitions: ", file_path)
	d = read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

	# Split the gene_id column (second) by the geneset_id column (first)
	gs = split(d[,2], d[,1])

	# Create the names list which is identity
	gs.names = as.list(names(gs))
	names(gs.names) = names(gs)

	# Create new GeneSet object
	object = new('GeneSet')

	# Populate the GeneSet object
	object@type = 'user-supplied'
	object@organism = 'user-supplied'
	object@dburl = 'user-supplied'

	object@set.gene = as.environment(gs)
	object@all.genes = as.character(unlist(gs, use.names = FALSE))

	object@set.name = as.environment(gs.names)
	message('Done setting up user-specified geneset..')

	return(object)
}

#' Function to read custom mappability files
#'
#' This function reads a two-columned tab-delimited text file (with header). Expected
#' column names are 'mappa' and 'gene_id'. Each line is for a unique 'gene_id'
#' and contains the mappability (between 0 and 1) for that gene.
#'
#' @param file_path A file path for the custom mappability.
#'
#' @return A \code{data.frame} containing \code{gene_id} and \code{mappa} columns.
read_mappa = function(file_path) {

	d = read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

	# Check columns.
	if(!(all(colnames(d) %in% c('mappa','gene_id')))) {
		stop("Error: header must contain columns named 'mappa' and 'gene_id'")
	}

	# Genes in this file should not be duplicated.
	if (sum(duplicated(d$gene_id)) > 0) {
		stop("Error reading mappability data: duplicate gene_ids exist in file.")
	}

	# Mappability should be between 0 and 1.
	if (min(d$mappa) < 0) {
		stop("Error: mappability must be >= 0.")
	}

	if (max(d$mappa) > 1) {
		stop("Error: mappability must be <= 1.")
	}

	return(d)
}
