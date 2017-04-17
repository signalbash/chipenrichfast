# Used in ..plot_polyenrich_spline(...)
avg_binned_numpeaks = function(gpw) {
	bygroup = stats::aggregate(cbind(num_peaks, length) ~ group, gpw, mean)
	bygroup$log_avg_length = log10(bygroup$length)
	names(bygroup) = c("group", "num_peaks", "avg_length", "log_avg_length")

	return(bygroup)
}

# Used in ..plot_polyenrich_spline(...)
calc_weights_polyenrich = function(gpw) {
	# Create model.
	model = "num_peaks ~ s(log10_length, bs = 'cr')"

	# Compute binomial spline fit.
	fit = gam(as.formula(model), data = gpw, family = "nb")

	fitted_numpeaks = fitted(fit)

	# These are not used, likely part of development
	# Perhaps of use for Chris in weighted polyenrich
	w0 = 1 / (fitted_numpeaks / mean(gpw$num_peaks, na.rm = TRUE))
	w0 = w0 / mean(w0, na.rm = TRUE)
	gpw$weight = w0
	gpw$resid.dev = resid(fit, type="deviance")

	# This is the only part that gets used
	gpw$fitted_numpeaks = fitted_numpeaks

	return(gpw)
}

#' QC plot for Poly-Enrich
#'
#' Create a plot the relationship between number of peaks assigned to a gene and
#' locus length. The plot shows an empirical fit to the data using a binomial
#' smoothing spline.
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
#' @param locusdef One of: 'nearest_tss', 'nearest_gene', 'exon', 'intron', '1kb',
#' '1kb_outside', '1kb_outside_upstream', '5kb', '5kb_outside', '5kb_outside_upstream',
#' '10kb', '10kb_outside', '10kb_outside_upstream'. For a description of each,
#' see the vignette or \code{\link{supported_locusdefs}}. Alternately, a file path for
#' a custom locus definition. NOTE: Must be for a \code{supported_genome()}, and
#' must have columns 'chr', 'start', 'end', and 'gene_id' or 'geneid'. For an
#' example custom locus definition file, see the vignette.
#' @param genome One of the \code{supported_genomes()}.
#' @param mappability One of \code{NULL}, a file path to a custom mappability file,
#' or an \code{integer} for a valid read length given by \code{supported_read_lengths}.
#' If a file, it should contain a header with two column named 'gene_id' and 'mappa'.
#' Gene IDs should be Entrez IDs, and mappability values should range from 0 and 1.
#' For an example custom mappability file, see the vignette. Default value is NULL.
#' @param legend If true, a legend will be drawn on the plot.
#' @param xlim Set the x-axis limit. NULL means select x-lim automatically.
#' @param ylim Set the y-axis limit. NULL means select y-lim automatically.
#'
#' @return A trellis plot object.
#'
#' @examples
#'
#' # Spline plot for E2F4 example peak dataset.
#' data(peaks_E2F4, package = 'chipenrich.data')
#'
#' # Create the plot for a different locus definition
#' # to compare the effect.
#' plot_polyenrich_spline(peaks_E2F4, locusdef = 'nearest_gene', genome = 'hg19')
#'
#' @export
#' @include constants.R utils.R supported.R setup.R randomize.R
#' @include read.R assign_peaks.R peaks_per_gene.R
plot_polyenrich_spline = function(peaks, locusdef = "nearest_tss", genome = supported_genomes(), mappability = NULL, legend = TRUE, xlim = NULL, ylim = NULL) {
	genome = match.arg(genome)

	ldef_list = setup_locusdef(locusdef, genome)
	ldef = ldef_list[['ldef']]
	tss = ldef_list[['tss']]

	mappa = setup_mappa(mappa_code = mappability, genome = genome, ldef_code = locusdef, ldef_obj = ldef)

	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks)
	} else if (class(peaks) == "character") {
		peakobj = read_bed(peaks)
	}
	num_peaks = length(peakobj)

	assigned_peaks = assign_peaks(peakobj, ldef, tss)
	ppg = num_peaks_per_gene(assigned_peaks, ldef, mappa)

	# Make plot.
	plotobj = ..plot_polyenrich_spline(gpw = ppg, mappability = mappability, num_peaks = num_peaks, legend = legend, xlim = xlim, ylim = ylim)
	return(plotobj)
}

..plot_polyenrich_spline = function(gpw, mappability, num_peaks, legend = TRUE, xlim = NULL, ylim = NULL) {
	############################################################################
	# Prepare the gpw table
	gpw = gpw[order(gpw$log10_length), ]
	gpw$group = ceiling((1:nrow(gpw)) / 25)

	############################################################################
	# Quantities for the plot

	# Scatterplot stuff
	avg_bins = avg_binned_numpeaks(gpw)

	# Spline stuff
	gpw = calc_weights_polyenrich(gpw) # gpw = genes, peaks, weights
	# genome_length = sum(as.numeric(gpw$length))
	# gpw$false_prob = 1 - (1 - (gpw$length / genome_length))^(num_peaks)

	############################################################################
	# Plotting parameters
	col_spline = "darkorange"

	panel_func = function(x,y,...) {
		lattice::panel.xyplot(x, y, ...)
	}

	custom_key = list(
		text = list(
		c(
			"Binomial Smoothing Spline Fit",
			"Avg Number of Peaks per Gene in Bins"
			)
		),
		lines = list(
			pch = c(17, 20),
			type = c("l", "p"),
			col = c(col_spline, "black"),
			lwd = 3
		),
		cex = 1.2
	)

	if (!legend) {
		custom_key = NULL
	}

	if (is.null(mappability)) {
		xlab = expression(paste(Log[10], " Locus Length"))
	} else {
		xlab = expression(paste(Log[10], " Mappable Locus Length"))
	}

	xmin_nopad = base::ifelse(is.null(xlim[1]), floor(min(gpw$log10_length)), floor(xlim[1]))
	xmax_nopad = base::ifelse(is.null(xlim[2]), ceiling(max(gpw$log10_length)), ceiling(xlim[2]))

	ymin_nopad = base::ifelse(is.null(ylim[1]), floor(min(avg_bins$num_peaks)), floor(ylim[1]))
	ymax_nopad = base::ifelse(is.null(ylim[2]), ceiling(max(avg_bins$num_peaks)), ceiling(ylim[2]))

	scales = list(
		x = list(
			axs = 'i',
			at = seq(xmin_nopad, xmax_nopad, 1)
		),
		y = list(
			axs = 'i',
			at = seq(ymin_nopad, ymax_nopad, 0.2)
		)
	)

	############################################################################
	# Plot
	plotobj = lattice::xyplot(
		fitted_numpeaks ~ log10_length,
		gpw,
		xlab = list(label = xlab, cex = 1.4),
		ylab = list(label = "Avg Number of Peaks per Gene in Bins", cex = 1.4),
		ylim = c(ymin_nopad - 0.5, ymax_nopad + 0.5),
		xlim = c(xmin_nopad - 0.5, xmax_nopad + 0.5),
		panel = panel_func,
		type = "l",
		lwd = 3,
		key = custom_key,
		scales = list(cex = 1.4),
		par.settings = lattice::simpleTheme(pch = c(17), col = c(col_spline))
	) + latticeExtra::as.layer(lattice::xyplot(num_peaks ~ log_avg_length, avg_bins, pch = 20, cex = 0.4, col = "black"))

	return(plotobj)
}
