# Used in ..plot_gene_coverage(...)
avg_binned_coverage = function(gpw) {
	bygroup = stats::aggregate(cbind(ratio, length) ~ group, gpw, mean)
	bygroup$log_avg_length = log10(bygroup$length)
	names(bygroup) = c("group", "ratio", "length", "log_avg_length")

	return(bygroup)
}

#' Plot QC plot for Broad-Enrich
#'
#' Create a plot showing the relationship between log10 locus length and the
#' proportion of gene loci covered by peaks.
#'
#' @param peaks Either a file path or a \code{data.frame} of peaks in BED-like
#' format. If a file path, the following formats are fully supported via their
#' file extensions: .bed, .broadPeak, .narrowPeak, .gff3, .gff2, .gff, and .bedGraph
#' or .bdg. BED3 through BED6 files are supported under the .bed extension. Files
#' without these extensions are supported under the conditions that the first 3
#' columns correspond to 'chr', 'start', and 'end' and that there is either no
#' header column, or it is commented out. If a \code{data.frame} A BEDX+Y style
#' \code{data.frame}. See \code{GenomicRanges::makeGRangesFromDataFrame} for
#' acceptable column names.
#' @param locusdef One of 'nearest_tss', 'nearest_gene', 'exon', 'intron', '1kb',
#' '1kb_outside', '1kb_outside_upstream', '5kb', '5kb_outside', '5kb_outside_upstream',
#' '10kb', '10kb_outside', '10kb_outside_upstream'. Alternately, a file path for
#' a custom locus definition. NOTE: Must be for a \code{supported_genome()}, and
#' must have columns 'chr', 'start', 'end', and 'gene_id' or 'geneid'.
#' @param genome One of the \code{supported_genomes()}.
#' @param mappability One of \code{NULL}, a file path to a custom mappability file,
#' or an \code{integer} for a valid read length given by \code{supported_read_lengths}.
#' If a file, it should contain a header with two column named 'gene_id' and 'mappa'.
#' Gene IDs should be Entrez IDs, and mappability values should range from 0 and 1.
#' Default value is NULL.
#' @param legend If true, a legend will be drawn on the plot.
#' @param xlim Set the x-axis limit. NULL means select x-lim automatically.
#'
#' @return A trellis plot object.
#'
#' @examples
#'
#' # Spline plot for E2F4 example peak dataset.
#' data(peaks_H3K4me3_GM12878, package = 'chipenrich.data')
#'
#' # Create the plot for a different locus definition
#' # to compare the effect.
#' plot_gene_coverage(peaks_H3K4me3_GM12878, locusdef = 'nearest_gene', genome = 'hg19')
#'
#' @export
#' @include constants.R utils.R supported.R setup.R randomize.R
#' @include read.R assign_peaks.R peaks_per_gene.R
plot_gene_coverage = function(peaks, locusdef = "nearest_tss", genome = supported_genomes(), mappability = NULL, legend = TRUE, xlim = NULL) {
	genome = match.arg(genome)

	ldef_list = setup_locusdef(locusdef, genome)
	ldef = ldef_list[['ldef']]

	mappa = setup_mappa(mappa_code = mappability, genome = genome, ldef_code = locusdef, ldef_obj = ldef)

	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks, genome = genome)
	} else if (class(peaks) == "character") {
		peakobj = read_bed(peaks, genome = genome)
	}
	num_peaks = length(peakobj)

	assigned_peaks = assign_peak_segments(peakobj, ldef)
	ppg = num_peaks_per_gene(assigned_peaks, ldef, mappa)
	ppg = calc_peak_gene_overlap(assigned_peaks, ppg)

	# Make plot.
	plotobj = ..plot_gene_coverage(gpw = ppg, mappability = mappability, num_peaks = num_peaks, xlim = xlim)
	return(plotobj)
}

..plot_gene_coverage = function(gpw, mappability, num_peaks, xlim = NULL) {
	############################################################################
	# Prepare the gpw table
	gpw = gpw[order(gpw$log10_length), ]
	gpw$group = ceiling((1:nrow(gpw)) / 25)

	############################################################################
	# Quantities for the plot

	# Scatterplot stuff
	avg_bins = avg_binned_coverage(gpw)

	############################################################################
	# Plotting parameters
	if (is.null(mappability)) {
		xlab = expression(paste(Log[10], " Locus Length"))
	} else {
		xlab = expression(paste(Log[10], " Mappable Locus Length"))
	}

	xmin_nopad = base::ifelse(is.null(xlim[1]), floor(min(gpw$log10_length)), floor(xlim[1]))
	xmax_nopad = base::ifelse(is.null(xlim[2]), ceiling(max(gpw$log10_length)), ceiling(xlim[2]))

	scales = list(
		x = list(
			axs = 'i',
			at = seq(xmin_nopad, xmax_nopad, 1)
		),
		y = list(
			axs = 'i',
			at = seq(0, 1, 0.2)
		)
	)

	############################################################################
	# Plot
	plotobj = lattice::xyplot(
		ratio ~ log_avg_length,
		avg_bins,
		main = 'Binned Locus Length versus Peak Coverage',
		xlab = list(label = xlab, cex = 1.4),
		ylab = list(label = "Proportion of Locus Covered by Peaks", cex = 1.4),
		ylim = c(-0.05, 1.05),
		xlim = c(xmin_nopad - 0.5, xmax_nopad + 0.5),
		pch = 20,
		cex = 0.4,
		col = 'black',
		scales = list(cex = 1.4)
	)

	return(plotobj)
}
