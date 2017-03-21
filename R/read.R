#' Read BEDX+Y files
#'
#' Given a \code{file_path}, read in BEDX+Y files. It is expected that files have
#' extensions that can reliably be parsed so that the correct variant of
#' \code{rtracklayer::import()} can be called. See \code{file_path} below.
#'
#' A warning is given if input regions are overlapping. In the case of enrichment
#' testing with \code{method = 'broadenrich'}, regions should be disjoint.
#'
#' Header rows must be commented with \code{#} to be ignored. Otherwise, and error
#' may result.
#'
#' Typically, this function will not be used alone, but inside \code{chipenrich()}.
#'
#' @param file_path A path to a file having the extensions: .bed, .broadPeak,
#' .narrowPeak, .gff3, .gff2, .gff, .wig, or .bedGraph or .bdg. Files may be
#' compressed with gzip, and so might end in .narrowPeak.gz.
#' @param genome The genome build used to get input regions. Usually inferred from
#' \code{chipenrich()}. Default is NA, in which case no genome is assigned to the
#' \code{GRanges} and the peaks cannot be sorted.
#'
#' @return A \code{GRanges} with \code{mcols} matching extra columns as needed.
#'
#' @examples
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
	# Establish extra columns and format as needed
	if(grepl('.narrowPeak', file_path)) {
		extraCols = c('signalValue' = 'numeric', 'pValue' = 'numeric', 'qValue' = 'numeric', 'peak' = 'integer')
		format = 'bed'
	} else if (grepl('.broadPeak', file_path)) {
		extraCols = c('signalValue' = 'numeric', 'pValue' = 'numeric', 'qValue' = 'numeric')
		format = 'bed'
	} else if (grepl('.gff3', file_path)) {
		extraCols = character()
		format = 'gff3'
	} else if (grepl('.gff2', file_path)) {
		extraCols = character()
		format = 'gff2'
	} else if (grepl('.gff', file_path)) {
		extraCols = character()
		format = 'gff'
	} else if (grepl('.wig', file_path)) {
		extraCols = character()
		format = 'wig'
	} else if (grepl('.bedGraph', file_path) || grepl('.bdg', file_path)) {
		extraCols = character()
		format = 'bedGraph'
	} else {
		extraCols = character()
		format = 'bed'
	}

	if(grepl('gff', format)) {
		# format = 'gff' doesn't use the extraCols parameter
		gr = rtracklayer::import(con = file_path, format = format, genome = genome)
	} else {
		gr = rtracklayer::import(con = file_path, format = format, genome = genome, extraCols = extraCols)
	}

	# Check for overlapping peaks
	reduced_gr = GenomicRanges::reduce(gr)
	if(length(gr) != length(reduced_gr)) {
		warning('Some input regions overlap. It is recommended that input regions be disjoint.')
	}

	# Sort
	gr = GenomicRanges::sort(gr)

	# Create name column, or replace what's there
	mcols(gr)$name = paste0('peak:', seq_along(gr))

	return(gr)
}

#' Convert a BEDX+Y data.frame and into GRanges
#'
#' Given a \code{data.frame} in BEDX+Y format, use the built-in function
#' \code{GenomicRanges::makeGRangesFromDataFrame()} to convert to \code{GRanges}.
#'
#' Typically, this function will not be used alone, but inside \code{chipenrich()}.
#'
#' @param dframe A BEDX+Y style \code{data.frame}.
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

	# Check for overlapping peaks
	reduced_gr = GenomicRanges::reduce(gr)
	if(length(gr) != length(reduced_gr)) {
		warning('Some input regions overlap. It is recommended that input regions be disjoint.')
	}

	# Sort
	gr = GenomicRanges::sort(gr)

	# Create name column, or replace what's there
	mcols(gr)$name = paste0('peak:', seq_along(gr))

	return(gr)
}
