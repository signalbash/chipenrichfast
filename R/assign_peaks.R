#' Assign peak midpoints to defined gene loci.
#'
#' Determine the midpoints of a set of input regions \code{peaks} and the overlap of the midpoints with a given locus definition \code{locusdef}. Also report the TSS that is nearest each region (peak) overlapping a defined locus and its distance.
#'
#' Typically, this function will not be used alone, but inside \code{chipenrich()}.
#'
#' @param peaks A \code{GRanges} object representing regions to be used for enrichment.
#' @param locusdef A locus definition object from \code{chipenrich.data}.
#' @param tss A \code{GRanges} object representing the TSSs for the genome build. Includes \code{mcols} for Entrez Gene ID \code{gene_id} and gene symbol \code{symbol}.
#'
#' @return A \code{data.frame} with columns for \code{peak_id, chr, peak_start, peak_end, gene_locus_start, gene_locus_end, gene_id, nearest_tss, nearest_tss_gene, dist_to_tss, nearest_tss_gene_strand}. The result is used in \code{num_peaks_per_gene()}.
#'
#' @examples
#'
#' data('locusdef.hg19.nearest_tss', package = 'chipenrich.data')
#' data('tss.hg19', package = 'chipenrich.data')
#'
#' file = system.file('extdata', 'test_assign.bed', package = 'chipenrich')
#' peaks = read_bed(file)
#'
#' assigned_peaks = assign_peaks(
#' 	peaks = peaks,
#' 	locusdef = locusdef.hg19.nearest_tss,
#' 	tss = tss.hg19)
#'
#' @export
assign_peaks = function(peaks, locusdef, tss, weighting=NULL) {
	# Extract GRanges of locusdef
	# NOTE: locusdef is an environment, so uses @
	ldef_gr = locusdef@granges

	# Determine midpoints of peaks and construct a GRanges object
	# on that basis. Include the peak name for later merging
	peak_mids = IRanges::mid(GenomicRanges::ranges(peaks))
	mids_gr = GenomicRanges::GRanges(
		seqnames = GenomeInfoDb::seqnames(peaks),
		ranges = IRanges::IRanges(start = peak_mids, end = peak_mids),
		name = GenomicRanges::mcols(peaks)$name
	)
    if (any(c("signalValue", "logsignalValue") %in% weighting)) {
        if (!("signalValue" %in% colnames(GenomicRanges::mcols(peaks)))) {
            stop("No signalValue column!")
        }
        mids_gr$signalValue = GenomicRanges::mcols(peaks)$signalValue
    }


	# Determine overlaps of peak midpoints with locus definition
	mid_ldef_overlaps = GenomicRanges::findOverlaps(mids_gr, ldef_gr)

	if(length(mid_ldef_overlaps) == 0) {
		stop('Intersection between peak midpoints and locus definition is empty. Try selecting a different locus definition.')
	}

	# Determine distance of peak midpoints to nearest TSSs
	mid_dist_to_ntss = GenomicRanges::distanceToNearest(mids_gr, tss)

	# Sign the distances by up (-) or down (+) stream
	mid_starts = GenomicRanges::start(mids_gr)[S4Vectors::queryHits(mid_dist_to_ntss)]
	tss_starts = GenomicRanges::start(tss)[S4Vectors::subjectHits(mid_dist_to_ntss)]
	tss_pos_strand = GenomicRanges::strand(tss)[S4Vectors::subjectHits(mid_dist_to_ntss)] == "+"
	tss_neg_strand = GenomicRanges::strand(tss)[S4Vectors::subjectHits(mid_dist_to_ntss)] == "-"

	neg1 = (mid_starts < tss_starts) & (tss_pos_strand)
	neg2 = (mid_starts > tss_starts) & (tss_neg_strand)
	neg_dists = as.logical(neg1 | neg2)

	# Reassign the signed distances
	dists = GenomicRanges::mcols(mid_dist_to_ntss)$distance
	dists[neg_dists] = dists[neg_dists] * -1
	GenomicRanges::mcols(mid_dist_to_ntss)$distance = dists

	# Construct data.frames for overlaps and distances
	mid_indices = S4Vectors::queryHits(mid_ldef_overlaps)
	ldef_indices = S4Vectors::subjectHits(mid_ldef_overlaps)
	mid_ldef_df = data.frame(
		peak_id = GenomicRanges::mcols(mids_gr)$name[mid_indices],
		chr = GenomeInfoDb::seqnames(mids_gr)[mid_indices],
		peak_start = GenomicRanges::start(peaks)[mid_indices],
		peak_end = GenomicRanges::end(peaks)[mid_indices],
		gene_locus_start = GenomicRanges::start(ldef_gr)[ldef_indices],
		gene_locus_end = GenomicRanges::end(ldef_gr)[ldef_indices],
		gene_id = GenomicRanges::mcols(ldef_gr)$gene_id[ldef_indices],
		gene_symbol = GenomicRanges::mcols(ldef_gr)$symbol[ldef_indices],
		stringsAsFactors = FALSE
	)
    
    # Add the signalValue column if method == "polyenrich_weighted"
    # Stops if there isn't one
    if (any(c("signalValue","logsignalValue") %in% weighting)) {
        mid_ldef_df$signalValue = GenomicRanges::mcols(peaks)$signalValue[mid_indices]
    }

	mid_indices = S4Vectors::queryHits(mid_dist_to_ntss)
	tss_indices = S4Vectors::subjectHits(mid_dist_to_ntss)
	mid_dist_df = data.frame(
		peak_id = GenomicRanges::mcols(mids_gr)$name[mid_indices],
		nearest_tss = GenomicRanges::start(tss)[tss_indices],
		nearest_tss_gene_id = GenomicRanges::mcols(tss)$gene_id[tss_indices],
		nearest_tss_symbol = GenomicRanges::mcols(tss)$symbol[tss_indices],
		dist_to_tss = GenomicRanges::mcols(mid_dist_to_ntss)$distance,
		nearest_tss_gene_strand = GenomicRanges::strand(tss)[tss_indices],
		stringsAsFactors = FALSE
	)

	# Merge the overlap and distance data.frames on the peak names
	d = base::merge(
		x = mid_ldef_df,
		y = mid_dist_df,
		by = 'peak_id',
		all.x = TRUE,
		all.y = FALSE,
		sort = FALSE
	)

	# Reorder the columns
	column_order = c(
		"peak_id",
		"chr",
		"peak_start",
		"peak_end",
		"gene_id",
		"gene_symbol",
		"gene_locus_start",
		"gene_locus_end",
		"nearest_tss",
		"dist_to_tss",
		"nearest_tss_gene_id",
		"nearest_tss_symbol",
		"nearest_tss_gene_strand"
	)
    if (any(c("signalValue","logsignalValue") %in% weighting)) {
        column_order = c(column_order, "signalValue")
    }
	d = d[, column_order]

	return(d)
}

#' Assign whole peaks to all overlapping defined gene loci.
#'
#' Determine all overlaps between the set of input regions \code{peaks} and the given locus definition \code{locusdef}. In addition, report where each overlap begins and ends, as well as the length of the overlap.
#'
#' Typically, this function will not be used alone, but inside \code{chipenrich()} with \code{method = 'broadenrich'}.
#'
#' @param peaks A \code{GRanges} object representing regions to be used for enrichment.
#' @param locusdef A locus definition object from \code{chipenrich.data}.
#'
#' @return A \code{data.frame} with columns for \code{peak_id, chr, peak_start, peak_end, gene_locus_start, gene_locus_end, gene_id, overlap_start, overlap_end, peak_overlap}. The result is used in \code{num_peaks_per_gene()}.
#'
#' @examples
#'
#' data('locusdef.hg19.nearest_tss', package = 'chipenrich.data')
#' data('tss.hg19', package = 'chipenrich.data')
#'
#' file = system.file('extdata', 'test_assign.bed', package = 'chipenrich')
#' peaks = read_bed(file)
#'
#' assigned_peaks = assign_peak_segments(
#' 	peaks = peaks,
#' 	locusdef = locusdef.hg19.nearest_tss)
#'
#' @export
assign_peak_segments = function(peaks, locusdef) {
	# Extract GRanges of locusdef
	# NOTE: locusdef is a class, so uses @
	ldef_gr = locusdef@granges

	# Determine overlaps of peaks with locus definition
	peak_ldef_overlaps = GenomicRanges::findOverlaps(peaks, ldef_gr)

	if(length(peak_ldef_overlaps) == 0) {
		stop('Intersection between peaks and locus definition is empty. Try selecting a different locus definition.')
	}

	# Determine start and ends for peaks and genes that overlap
	peak_indices = S4Vectors::queryHits(peak_ldef_overlaps)
	ldef_indices = S4Vectors::subjectHits(peak_ldef_overlaps)

	peak_start = GenomicRanges::start(peaks)[peak_indices]
	peak_end = GenomicRanges::end(peaks)[peak_indices]
	gene_start = GenomicRanges::start(ldef_gr)[ldef_indices]
	gene_end = GenomicRanges::end(ldef_gr)[ldef_indices]

	# Use pmin() and pmax() to determine componentwise mins and maxes
	# which are also the starts and ends of the overlaps
	overlap_start = pmax(peak_start, gene_start)
	overlap_end = pmin(peak_end, gene_end)
	overlap_length = overlap_end - overlap_start + 1

	d = data.frame(
		peak_id = GenomicRanges::mcols(peaks)$name[peak_indices],
		chr = GenomeInfoDb::seqnames(peaks)[peak_indices],
		peak_start = peak_start,
		peak_end = peak_end,
		gene_locus_start = gene_start,
		gene_locus_end = gene_end,
		gene_id = GenomicRanges::mcols(ldef_gr)$gene_id[ldef_indices],
		gene_symbol = GenomicRanges::mcols(ldef_gr)$symbol[ldef_indices],
		overlap_start = overlap_start,
		overlap_end = overlap_end,
		peak_overlap = overlap_length,
		stringsAsFactors = FALSE
	)

	return(d)
}
