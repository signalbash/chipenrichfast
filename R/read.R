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
#' for example. Forfiles with extension support, the \code{rtracklayer::import()}
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
