# Used in ..plot_gene_coverage(...)
avg_binned_coverage = function(gpw, bin_size = 25) {
	d = gpw
	d = d[order(d$log10_length), ]
	d$group = ceiling((1:dim(d)[1])/bin_size)

	bygroup = stats::aggregate(cbind(ratio, log10_length) ~ group, d, mean)
	names(bygroup) = c("group", "ratio", "log_avg_length")
	bygroup
}

#' Plot probability of peak being assigned to a gene vs. gene length
#'
#' Create a plot showing the probability of a gene being assigned a peak given
#' its locus length. The plot shows an empirical fit to the data using a binomial
#' smoothing spline.
#'
#' @param peaks A \code{data.frame}, or tab-delimited text file (BED, narrowPeak,
#' broadPeak, etc) with the first three columns being chrom, start, and end. The
#' data frame should have at least 3 columns: chrom, start, and end. Chrom
#' should follow UCSC convention, e.g. "chrX".
#' @param locusdef A string denoting the gene locus definition to be used, or
#' the full path to a user-defined locus definition file. A gene locus definition
#' controls how peaks are assigned to genes. See \code{\link{supported_locusdefs}}
#' for a list of supported definitions built-in. If using a user-specified file,
#' the file must have 4 columns: gene_id, chrom, start, end and be tab-delimited.
#' @param genome A string indicating the genome upon which the peaks file is
#' based. Supported genomes are listed by the \code{\link{supported_genomes}} function.
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
#' data(peaks_E2F4, package = 'chipenrich.data')
#' peaks_E2F4 = subset(peaks_E2F4, peaks_E2F4$chrom == 'chr1')
#' plot_gene_coverage(peaks_E2F4, genome = 'hg19')
#'
#' @export
#' @include constants.R utils.R supported.R setup.R randomize.R
#' @include read.R assign_peaks.R peaks_per_gene.R
plot_gene_coverage = function(peaks, locusdef = "nearest_tss", genome = supported_genomes(), mappability = NULL, legend = TRUE, xlim = NULL) {
	genome = match.arg(genome)

	ldef_list = setup_locusdef(locusdef, genome)
	ldef = ldef_list[['ldef']]

	mappa = setup_mappa(mappa_code = mappability, genome = genome, ldef_code = locusdef, ldef_obj = ldef)

	# Get peaks from user's file.
	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks, genome = genome)
	} else if (class(peaks) == "character") {
		peakobj = read_bed(peaks, genome = genome)
	}

	# Assign peaks to genes.
	assigned_peaks = assign_peak_segments(peakobj, ldef)
	peak_genes = unique(assigned_peaks$gene_id)

	ppg = num_peaks_per_gene(assigned_peaks, ldef, mappa=NULL)
	ppg = calc_peak_gene_overlap(assigned_peaks, ppg)

	# Make plot.
	plotobj = ..plot_gene_coverage(ppg)
	return(plotobj)
}

..plot_gene_coverage = function(ppg) {

	avg_bins = avg_binned_coverage(ppg, bin_size = 25)

	plotobj = xyplot(
		ratio ~ log_avg_length,
		avg_bins,
		main = 'Binned Locus Length versus Peak Coverage',
		xlab = 'log10(locus length)',
		ylab = 'Proportion of locus covered by peak',
		pch = 20,
		col = 'black'
	)

	return(plotobj)
}
