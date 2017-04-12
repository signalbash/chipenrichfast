# Used in ..plot_chipenrich_spline(...)
avg_binned_peak = function(gpw) {
	bygroup = stats::aggregate(cbind(peak, length) ~ group, gpw, mean)
	bygroup$log_avg_length = log10(bygroup$length)
	names(bygroup) = c("group", "peak", "avg_length", "log_avg_length")

	return(bygroup)
}

# Used in ..plot_chipenrich_spline(...)
calc_weights_chipenrich = function(gpw) {
	# Create model.
	model = "peak ~ s(log10_length, bs = 'cr')"

	# Compute binomial spline fit.
	fit = gam(as.formula(model), data = gpw, family = "binomial")

	# Compute weights for each gene, based on the predicted prob(peak) for each gene.
	ppeak = fitted(fit)

	# These are not used, likely part of development
	# Perhaps of use for Chris in weighted polyenrich
	w0 = 1 / (ppeak / mean(gpw$peak, na.rm = TRUE))
	w0 = w0 / mean(w0, na.rm = TRUE)
	gpw$weight = w0
	gpw$resid.dev = resid(fit, type="deviance")

	# This is the only part that gets used
	gpw$prob_peak = ppeak

	return(gpw)
}

#' QC plot for ChIP-Enrich
#'
#' A plot showing an approximation of the empirical spline used to model the
#' relationship between a gene having a peak and the locus length. For visual
#' clarity, genes are binned into groups of 25 after sorting by locus length.
#' Expected fits assuming independence of locus length and presence of a peak,
#' and assuming proportionality of locus length and presence of a peak are given
#' to demonstrate deviation from either for the dataset.
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
#' plot_chipenrich_spline(peaks_E2F4, locusdef = 'nearest_gene', genome = 'hg19')
#'
#' @export
#' @include constants.R utils.R supported.R setup.R randomize.R
#' @include read.R assign_peaks.R peaks_per_gene.R
plot_chipenrich_spline = function(peaks, locusdef = "nearest_tss", genome = supported_genomes(), mappability = NULL, legend = TRUE, xlim = NULL) {
	genome = match.arg(genome)

	ldef_list = setup_locusdef(locusdef, genome)
	ldef = ldef_list[['ldef']]
	tss = ldef_list[['tss']]

	mappa = setup_mappa(mappa_code = mappability, genome = genome, ldef_code = locusdef, ldef_obj = ldef)

	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks, genome = genome)
	} else if (class(peaks) == "character") {
		peakobj = read_bed(peaks, genome = genome)
	}
	num_peaks = length(peakobj)

	assigned_peaks = assign_peaks(peakobj, ldef, tss)
	ppg = num_peaks_per_gene(assigned_peaks, ldef, mappa)

	# Make plot.
	plotobj = ..plot_chipenrich_spline(gpw = ppg, mappability = mappability, num_peaks = num_peaks, legend = legend, xlim = xlim)
	return(plotobj)
}

..plot_chipenrich_spline = function(gpw, mappability, num_peaks, legend = TRUE, xlim = NULL) {
	############################################################################
	# Prepare the gpw table
	gpw = gpw[order(gpw$log10_length), ]
	gpw$group = ceiling((1:nrow(gpw)) / 25)

	############################################################################
	# Quantities for the plot

	# Scatterplot stuff
	avg_bins = avg_binned_peak(gpw)

	# Spline stuff
	gpw = calc_weights_chipenrich(gpw) # gpw = genes, peaks, weights
	genome_length = sum(as.numeric(gpw$length))
	gpw$false_prob = 1 - (1 - (gpw$length / genome_length))^(num_peaks)

	############################################################################
	# Plotting parameters
	col_rand_gene = "grey35"
	col_rand_peak = "grey74"
	col_spline = "darkorange"

	panel_func = function(x,y,...) {
		lattice::panel.xyplot(x, y, ...)
		lattice::panel.abline(h = mean(y), col = col_rand_gene, lwd = 3)
	}

	custom_key = list(
		text = list(
		c(
			"Expected Fit - Peaks Independent of Locus Length",
			"Expected Fit - Peaks Proportional to Locus Length",
			"Binomial Smoothing Spline Fit",
			"Proportion of Genes in Bin with at Least 1 Peak"
			)
		),
		lines = list(
			pch = c(20, 15, 17, 20),
			type = c("l", "l", "l", "p"),
			col = c(col_rand_gene, col_rand_peak, col_spline, "black"),
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
		false_prob + prob_peak ~ log10_length,
		gpw,
		xlab = list(label = xlab, cex = 1.4),
		ylab = list(label = "Prop. of Genes with at Least 1 Peak", cex = 1.4),
		ylim = c(-0.05, 1.05),
		xlim = c(xmin_nopad - 0.5, xmax_nopad + 0.5),
		panel = panel_func,
		type = "l",
		lwd = 3,
		key = custom_key,
		scales = list(cex = 1.4),
		par.settings = lattice::simpleTheme(pch = c(15,17), col = c(col_rand_peak, col_spline))
	) + latticeExtra::as.layer(lattice::xyplot(peak ~ log_avg_length, avg_bins, pch = 20, cex = 0.4, col = "black"))

	return(plotobj)
}
