#!/usr/bin/env Rscript
os = Sys.info()[1]

#' Run ChIP-Enrich on a dataset of ChIP-seq peaks
#'
#' Run gene set enrichment testing (ChIP-Enrich or Broad-Enrich) on a ChIP-seq
#' peak dataset or other type of dataset consisting of regions across the genome.
#' The user can call \code{chipenrich()} to run the method on their data. A number
#' of arguments can be provided to change the type of test, the genome build,
#' which sets of genes to test, how peaks are assigned to genes, and other minor options.
#'
#' @param peaks A \code{data.frame}, or tab-delimited text file (BED, narrowPeak,
#' broadPeak, etc) with the first three columns being chrom, start, and end. The
#' data frame should have at least 3 columns: chrom, start, and end. Chrom
#' should follow UCSC convention, e.g. "chrX".
#' @param out_name Prefix string to use for naming output files. This should not
#' contain any characters that would be illegal for the system being used (Unix,
#' Windows, etc.) The default value is "chipenrich", and a file "chipenrich_results.tab"
#' is produced. If \code{qc_plots} is set, then a file "chipenrich_qcplots.pdf"
#' is produced containing a number of quality control plots. If \code{out_name}
#' is set to NULL, no files are written, and results then must be retrieved from
#' the list returned by \code{chipenrich}.
#' @param out_path Directory to which results files will be written out. Defaults
#' to the current working directory as returned by \code{\link{getwd}}.
#' @param genome One of the \code{supported_genomes()}.
#' @param genesets A character vector of geneset databases to be tested for
#' enrichment. See \code{supported_genesets()}. Alternately, a file path to a
#' a tab-delimited text file with header and first column being the geneset ID
#' or name, and the second column being Entrez Gene IDs.
#' @param locusdef One of 'nearest_tss', 'nearest_gene', 'exon', 'intron', '1kb',
#' '1kb_outside', '1kb_outside_upstream', '5kb', '5kb_outside', '5kb_outside_upstream',
#' '10kb', '10kb_outside', '10kb_outside_upstream'. Alternately, a file path for
#' a custom locus definition. NOTE: Must be for a \code{supported_genome()}, and
#' must have columns 'chr', 'start', 'end', and 'gene_id', or 'geneid'.
#' @param method A character string specifying the method to use for enrichment
#' testing. Must be one of ChIP-Enrich ('chipenrich') (default), ChIP-Enrich fast
#' ('chipenrich_fast'), Broad-Enrich ('broadenrich'), Poly-Enrich ('polyenrich'),
#' Poly-Enrich fast ('polyenrich_fast'), or Fisher's exact test ('fet').
#' For a list of supported methods, use \code{\link{supported_methods}}.
#' @param fisher_alt If method is 'fet', this option indicates the alternative
#' for Fisher's exact test, and must be one of 'two-sided' (default), 'greater',
#' or 'less'.
#' @param mappability One of \code{NULL}, a file path to a custom mappability file,
#' or an \code{integer} for a valid read length given by \code{supported_read_lengths}.
#' If a file, it should contain a header with two column named 'gene_id' and 'mappa'.
#' Gene IDs should be Entrez IDs, and mappability values should range from 0 and 1.
#' Default value is NULL.
#' @param qc_plots A logical variable that enables the automatic generation of
#' plots for quality control.
#' @param min_geneset_size Sets the minimum number of genes a gene set may have
#' to be considered for enrichment testing.
#' @param max_geneset_size Sets the maximum number of genes a gene set may have
#' to be considered for enrichment testing.
#' @param num_peak_threshold Sets the threshold for how many peaks a gene must
#' have to be considered as having a peak. Defaults to 1. Only relevant for
#' Fisher's exact test and ChIP-Enrich methods.
#' @param n_cores The number of cores to use for enrichment testing. We recommend
#' using only up to the maximum number of \emph{physical} cores present, as
#' virtual cores do not significantly decrease runtime. Default number of cores
#' is set to 1. NOTE: Windows does not support multicore enrichment.
#'
#' @return A list, containing the following items:
#' \item{peaks }{
#' A data frame containing peak assignments to genes. Peaks which do not overlap
#' a gene locus are not included. Each peak that was assigned to a gene is listed,
#' along with the peak midpoint or peak interval coordinates (depending on which
#' was used), the gene to which the peak was assigned, the locus start and end
#' position of the gene, and the distance from the peak to the TSS.
#'
#' The columns are:
#'
#' \describe{
#'   \item{peak_id}{ is an ID given to unique combinations of chromosome, peak
#' start, and peak end. }
#'   \item{chrom}{ is the chromosome the peak originated from. }
#'   \item{peak_start}{ is start position of the peak. }
#'   \item{peak_end}{ is end position of the peak. }
#'   \item{peak_midpoint}{ is the midpoint of the peak. }
#'   \item{gene_id}{ is the Entrez ID of the gene to which the peak was assigned. }
#'   \item{gene_symbol}{ is the official gene symbol for the gene_id (above). }
#'   \item{gene_locus_start}{ is the start position of the locus for the gene to which the peak was assigned (specified by the locus definition used.) }
#'   \item{gene_locus_end}{ is the end position of the locus for the gene to which the peak was assigned (specified by the locus definition used.) }
#'   \item{nearest_tss}{ (\code{method='chipenrich'} and \code{method='fet'}) is the closest TSS to this peak (for any gene, not necessarily the gene this peak was assigned to.) }
#'   \item{nearest_tss_gene}{ (\code{method='chipenrich'} and \code{method='fet'}) is the gene having the closest TSS to the peak (should be the same as gene_id when using the nearest TSS locus definition.) }
#'   \item{nearest_tss_gene_strand}{ (\code{method='chipenrich'} and \code{method='fet'}) is the strand of the gene with the closest TSS. }
#'   \item{overlap_start}{ (\code{method='broadenrich'} only) the start position of the peak overlap with the gene locus.}
#'   \item{overlap_end}{ (\code{method='broadenrich'} only) the end position of the peak overlap with the gene locus.}
#'   \item{peak_overlap}{ (\code{method='broadenrich'} only) the base pair overlap of the peak with the gene locus.}
#' }}
#'
#' \item{results }{
#' A data frame of the results from performing the gene set enrichment test on
#' each geneset that was requested (all genesets are merged into one final data
#' frame.) The columns are:
#'
#' \describe{
#'   \item{Geneset.ID}{ is the identifier for a given gene set from the selected database.  For example, GO:0000003. }
#'   \item{Geneset.Type}{ specifies from which database the Geneset.ID originates.  For example, "Gene Ontology Biological Process."}
#'   \item{Description}{ gives a definition of the geneset. For example, "reproduction."}
#'   \item{P.Value}{ is the probability of observing the degree of enrichment of the gene set given the null hypothesis that peaks are not associated with any gene sets.}
#'   \item{FDR}{ is the false discovery rate proposed by Bejamini \& Hochberg for adjusting the p-value to control for family-wise error rate.}
#'   \item{Odds.Ratio}{ is the estimated odds that peaks are associated with a given gene set compared to the odds that peaks are associated with other gene sets, after controlling for locus length and/or mappability.  An odds ratio greater than 1 indicates enrichment, and less than 1 indicates depletion.}
#'   \item{N.Geneset.Genes}{ is the number of genes in the gene set.}
#'   \item{N.Geneset.Peak.Genes}{ is the number of genes in the genes set that were assigned at least one peak.}
#'   \item{Geneset.Avg.Gene.Length}{ is the average length of the genes in the gene set.}
#'   \item{Geneset.Avg.Gene.Coverage}{ (\code{method='broadenrich'} only) is the mean proportion of the gene loci in the gene set covered by a peak.}
#'   \item{Geneset.Peak.Genes}{ is the list of genes from the gene set that had at least one peak assigned.}
#'
#' }}
#'
#' \item{opts }{A data frame containing the arguments/values passed to \code{chipenrich}.}
#'
#' \item{peaks_per_gene }{
#' A data frame of the count of peaks per gene. The columns are:
#'
#' \describe{
#'   \item{gene_id}{ is the Entrez Gene ID. }
#'   \item{length}{ is the length of the gene's locus (depending on which locus definition you chose.)}
#'   \item{log10_length}{ is the log10(locus length) for the gene.}
#'   \item{num_peaks}{ is the number of peaks that were assigned to the gene, given the current locus definition. }
#'   \item{peak}{ is whether or not the gene is considered to have a peak, as defined by \code{num_peak_threshold}. }
#'   \item{peak_overlap}{ (\code{method='broadenrich'} only) is the number of base pairs of the gene covered by a peak.}
#'   \item{ratio}{ (\code{method='broadenrich'} only) is the proportion of the gene covered by a peak.}
#' }}
#'
#' @examples
#'
#' # Run ChipEnrich using an example dataset, assigning peaks to the nearest TSS,
#' # testing all Biocarta and Panther pathways
#' data(peaks_E2F4, package = 'chipenrich.data')
#' peaks_E2F4 = subset(peaks_E2F4, peaks_E2F4$chrom == 'chr1')
#' gs_path = system.file('extdata','vignette_genesets.txt', package='chipenrich')
#' results = chipenrich(peaks_E2F4, method='chipenrich', locusdef='nearest_tss',
#' 			genome = 'hg19', genesets=gs_path, out_name=NULL)
#'
#' # Get the list of peaks that were assigned to genes.
#' assigned_peaks = results$peaks
#'
#' # Get the results of enrichment testing.
#' enrich = results$results
#'
#' @export
#' @include constants.R utils.R supported.R setup.R randomize.R
#' @include read.R assign_peaks.R peaks_per_gene.R
#' @include plot_dist_to_tss.R plot_gene_coverage.R plot_spline_length.R
#' @include test_approx.R test_binomial.R test_fisher.R test_gam.R
chipenrich = function(
	peaks,
	out_name = "chipenrich",
	out_path = getwd(),
	genome = supported_genomes(),
	genesets = c(
		'GOBP',
		'GOCC',
		'GOMF'),
	locusdef = "nearest_tss",
	method = 'chipenrich',
	mappability = NULL,
	fisher_alt = "two.sided",
	qc_plots = TRUE,
	min_geneset_size = 15,
	max_geneset_size = 2000,
	num_peak_threshold = 1,
	n_cores = 1
) {
	genome = match.arg(genome)

	############################################################################
	# Collect options for opts output
	opts_list = as.list(sys.call())
	opts_list = opts_list[2:length(opts_list)]

	opts = data.frame(
		parameters = names(opts_list),
		values = as.character(opts_list),
		stringsAsFactors = FALSE
	)

	############################################################################
	# Deal with randomizations if present
	# Randomizations are accessed by appending _rndall or _rndlength to
	# the geneset names. Parsing is done here.
	rndall = all(grepl('rndall',genesets))
	rndlength = all(grepl('rndlength',genesets))
	rndloc = all(grepl('rndloc',genesets))

	if(rndall) {
		genesets = gsub('_rndall','',genesets)
	} else if (rndlength) {
		genesets = gsub('_rndlength','',genesets)
	} else if (rndloc) {
		genesets = gsub('_rndloc','',genesets)
	}

	############################################################################
	############################################################################
	# Warnings based on OS or method/ldef combinations
	############################################################################
	############################################################################

	############################################################################
	# CHECK OS and give multicore warning if Windows
	if(os == 'Windows') {
		# NOTE: instead of warning, set n_cores to 1 and remove all OS checks
		# in the test_*() functions.
		message('Warning! Multicore enrichment is not supported on Windows.')
	}

	############################################################################
	############################################################################
	# Checks and genome, locusdef, geneset, and mappa setup
	############################################################################
	############################################################################

	ldef_list = setup_locusdef(locusdef, genome, rndloc)
	ldef = ldef_list[['ldef']]
	tss = ldef_list[['tss']]
	if(rndloc) {
		ldef = randomize_locusdef(ldef, 50)
	}

	geneset_list = setup_genesets(gs_codes = genesets, ldef_obj = ldef, genome = genome, min_geneset_size = min_geneset_size, max_geneset_size = max_geneset_size)

	mappa = setup_mappa(mappa_code = mappability, genome = genome, ldef_code = locusdef, ldef_obj = ldef)

	############################################################################
	# CHECK method and get() it if okay
	get_test_method = function(x) {
		if (method %in% names(SUPPORTED_METHODS)) {
			return(SUPPORTED_METHODS[[method]])
		} else if (method %in% names(HIDDEN_METHODS)) {
			return(HIDDEN_METHODS[[method]])
		} else {
			stop(sprintf("Error: invalid enrichment test requested: %s, contact developer.",method))
		}
	}
	testf = get_test_method(method)
	test_func = get(testf)

	# Test name.
	method_name = METHOD_NAMES[[method]]

	############################################################################
	# Warn user if they are trying to use FET with a
	# locus definition that might lead to biased results.
	if (method == "fet") {
		if (is.character(locusdef) && !locusdef %in% c("1kb","5kb")) {
			message("Warning: Fisher's exact test should only be used with the 1kb or 5kb locus definition.")
		}
	}

	# Warn user if they are using the binomial test.
	if (method == "binomial") {
		message("Warning: the binomial test is provided for comparison purposes only.")
		message("This test will almost always give biased results favoring gene sets with short average locus length.")
	}

	############################################################################
	############################################################################
    # Start enrichment process
	############################################################################
	############################################################################

	######################################################
	# Read in and format peaks (from data.frame or file)
	if (class(peaks) == "data.frame") {
		message('Reading peaks from data.frame...')
		peakobj = load_peaks(peaks, genome = genome)
	} else if (class(peaks) == "character") {
		peakobj = read_bed(peaks, genome = genome)
	}

	# Number of peaks in data.
	num_peaks = length(peakobj)

	######################################################
	# Assign peaks to genes. NOTE: If method = 'broadenrich' use
	# assign_peak_segments(), otherwise use assign_peaks().
	if(!(method == 'broadenrich' || method == 'broadenrich_splineless')) {
		message("Assigning peaks to genes with assign_peaks(...) ..")
		assigned_peaks = assign_peaks(peakobj, ldef, tss)
		message("Successfully assigned peaks..")
	} else {
		message("Assigning peaks to genes with assigned_peak_segments(...) ..")
		assigned_peaks = assign_peak_segments(peakobj, ldef)
		message("Successfully assigned peaks..")
	}

	peak_genes = unique(assigned_peaks$gene_id)

	######################################################
	# Compute peaks per gene table
	ppg = num_peaks_per_gene(assigned_peaks, ldef, mappa)
	# This seems redundant given num_peaks_per_gene(...)
	ppg$peak = recode_peaks(ppg$num_peaks, num_peak_threshold)

	# Add relevant columns to ppg depending on the method
	if(method == 'broadenrich' || method == 'broadenrich_splineless') {
		message("Calculating peak overlaps with gene loci..")
		ppg = calc_peak_gene_overlap(assigned_peaks,ppg)
	}
	if(method == 'chipapprox') {
		message("Calculating weights for approximate method..")
		ppg = calc_approx_weights(ppg,mappa)
	}

	######################################################
	# Randomize ppg table if randomization API invoked
	# Catch randomizations if present
	if(rndall) {
		message('Randomizing across all genes.')
		ppg = randomize_ppg_all(ppg)
	} else if (rndlength) {
		message('Randomizing within length bins.')
		ppg = randomize_ppg_length(ppg)
	}

	######################################################
	# Enrichment
	# Run chipenrich method on each geneset.
	results = list()
	for (gobj in geneset_list) {
		message(sprintf("Test: %s",method_name))
		message(sprintf("Genesets: %s",gobj@type))
		message("Running tests..")
		if (testf == "test_gam") {
			rtemp = test_func(gobj,ppg,n_cores)
		}
		if (testf == "test_fisher_exact") {
			rtemp = test_func(gobj,ppg,alternative=fisher_alt)
		}
		if (testf == "test_binomial") {
			rtemp = test_func(gobj,ppg)
		}
		if (testf == "test_gam_ratio") {
			rtemp = test_func(gobj,ppg,n_cores)
		}
		if (testf == "test_approx") {
			rtemp = test_func(gobj,ppg,nwp=FALSE,n_cores)
		}
		if (testf == "test_gam_ratio_splineless") {
			rtemp = test_func(gobj,ppg,n_cores)
		}
		if (testf == "test_gam_nb") {
			rtemp = test_func(gobj,ppg,n_cores);
		}
		if (testf == "test_gam_nb_fast") {
			rtemp = test_func(gobj,ppg,n_cores);
		}
		if (testf == "test_gam_fast") {
			rtemp = test_func(gobj,ppg,n_cores);
		}

		# Annotate with geneset descriptions.
		rtemp$"Description" = as.character(mget(rtemp$Geneset.ID, gobj@set.name, ifnotfound=NA))
		rtemp$"Geneset.Type" = gobj@type;

		results[[gobj@type]] = rtemp;
	}
	enrich = Reduce(rbind,results)

	######################################################
	# Post-process enrichment with various orderings

	# Re-order the columns to something sensible.
	column_order = c(
		"Geneset.Type",
		"Geneset.ID",
		"Description",
		"P.value",
		"FDR",
		"Effect",
		"Odds.Ratio",
		"P.Success",
		"Status",
		"N.Geneset.Genes",
		"N.Geneset.Peak.Genes",
		"Geneset.Avg.Gene.Length",
		"Geneset.Avg.Gene.Coverage",
		"Geneset.Peak.Genes")
	column_order = intersect(column_order, names(enrich))
	enrich = enrich[, column_order]

	# Order results by p-value.
	enrich = enrich[order(enrich$P.value), ]

	# If there is a status column, re-sort so enriched terms are on top.
	if ("Status" %in% names(enrich)) {
		enrich = enrich[order(enrich$Status, decreasing=TRUE), ]
	}

	# Pull out tests that failed.
	bad_enrich = subset(enrich, is.na(enrich$P.value))
	enrich = subset(enrich, !is.na(enrich$P.value))
	rownames(enrich) = c(1:nrow(enrich))

	######################################################
	# Write result objects to files
	if (!is.null(out_name)) {
		filename_analysis = file.path(out_path, sprintf("%s_results.tab", out_name))
		write.table(enrich, file = filename_analysis, row.names = FALSE, quote = FALSE, sep = "\t")
		message("Wrote results to: ", filename_analysis)

		filename_peaks = file.path(out_path, sprintf("%s_peaks.tab", out_name))
		write.table(assigned_peaks, file = filename_peaks, row.names = FALSE, quote = FALSE, sep = "\t")
		message("Wrote peak-to-gene assignments to: ", filename_peaks)

		filename_opts = file.path(out_path, sprintf("%s_opts.tab", out_name))
		write.table(opts, file = filename_opts, row.names = FALSE, quote = FALSE, sep = "\t")
		message("Wrote run options/arguments to: ", filename_opts)

		filename_ppg = file.path(out_path, sprintf("%s_peaks-per-gene.tab", out_name))
		write.table(ppg, file = filename_ppg, row.names = FALSE, quote = FALSE, sep = "\t")
		message("Wrote count of peaks per gene to: ", filename_ppg)

		# If the user requested QC plots, generate those as well.
		# -- Spline fit plot
		# -- Histogram of distance from peaks to TSSs
		# -- Histogram of p-values from test
		# -- Expected # peaks vs. observed # peaks
		if (qc_plots) {
			filename_qcplots = file.path(out_path, sprintf("%s_qcplots.pdf", out_name))
			grDevices::pdf(filename_qcplots)
				if (!(method=='broadenrich' || method=='broadenrich_splineless')) {
					print(..plot_spline_length(ldef, peak_genes, num_peaks, mappa=mappa))
					print(..plot_dist_to_tss(peakobj, tss))
				} else {
					print(..plot_gene_coverage(ppg))
				}
			grDevices::dev.off()
			message("Wrote QC plots to: ",filename_qcplots)
		}
	}

	######################################################
	# Return objects as list
	return(list(
		peaks = assigned_peaks,
		results = enrich,
		opts = opts,
		peaks_per_gene = ppg
	))
}
