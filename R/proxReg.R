#' Run Proximity Regulation test on narrow genomic regions
#' 
#' Description
#' 
#' @section Method:
#' ..
#' 
#' @section ...
#' 
#' @param xx blah
#' @param xx2 blah
#' ..
#' 
#' @return ...
#' 
#' @examples 
#' ...
#' 
#' @export
#' @include constants.R utils.R supported.R setup.R randomize.R
proxReg = function(
	peaks,
	out_name = "proxReg",
	out_path = getwd(),
	genome = supported_genomes(),
	reglocation = "tss",
	genesets = c(
		'GOBP',
		'GOCC',
		'GOMF'),
	qc_plots = TRUE,
	min_geneset_size = 15,
	max_geneset_size = 2000,
	n_cores = 1
) {
	genome = match.arg(genome)
	
	n_cores = reset_ncores_for_windows(n_cores)
	
	# Kill if reglocation is not tss or enhancer
	if (!(reglocation %in% c("tss","enhancer"))) {
		stop("Unsupported regulatory location!")
	}
	
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
	# Setup locus definitions, genesets, and mappability
	
	#Locus definition will always be NTSS for this
	ldef_list = setup_locusdef(ldef_code = "nearest_tss", genome, randomization = NULL)
	ldef = ldef_list[['ldef']]
	tss = ldef_list[['tss']]
	
	geneset_list = setup_genesets(gs_codes = genesets, ldef_obj = ldef, genome = genome, min_geneset_size = min_geneset_size, max_geneset_size = max_geneset_size)
	
	############################################################################
	############################################################################
	# Start enrichment process
	############################################################################
	############################################################################
	
	######################################################
	# Read in and format peaks (from data.frame or file)
	if (class(peaks) == "data.frame") {
		message('Reading peaks from data.frame...')
		peakobj = load_peaks(peaks)
	} else if (class(peaks) == "character") {
		peakobj = read_bed(peaks)
	}
	
	# Number of peaks in data.
	num_peaks = length(peakobj)
	
	######################################################
	# Assign peaks to genes.
	message("Assigning peaks to genes with assign_peaks(...) ..")
	assigned_peaks = assign_peaks(peakobj, ldef, tss)
	
	# Creating a column for adjusted DTSS
	if ("tss" %in% reglocation) {
		assigned_peaks$log_dtss = log(abs(assigned_peaks$dist_to_tss)+1)
		log_gene_ll = log(assigned_peaks$gene_locus_end-assigned_peaks$gene_locus_start)
		pred_log_dtss = as.numeric(fitted(chipenrich.data::spline.log_dtss.90ENCODE, peaks, type="terms"))
		assigned_peaks$scaled_dtss = assigned_peaks$log_dtss-pred_log_dtss
	}
	
	# Creating a column for enhancer distances
	if ("enhancer" %in% reglocation) {
		enhancers = chipenrich.data::enhancer.dnase_thurman.0
		peakobj2 = GenomicRanges::makeGRangesFromDataFrame(assigned_peaks,
					seqnames.field = "chr", start.field = "peak_start",end.field = "peak_end")
		peak_mids = IRanges::mid(GenomicRanges::ranges(peakobj2))
		mids_gr = GenomicRanges::GRanges(
			seqnames = GenomeInfoDb::seqnames(peakobj2),
			ranges = IRanges::IRanges(start = peak_mids, end = peak_mids),
			name = GenomicRanges::mcols(peakobj2)$name
		)
		enhancer_mids = IRanges::mid(GenomicRanges::ranges(enhancers))
		enhancer_mids_gr = GenomicRanges::GRanges(
			seqnames = GenomeInfoDb::seqnames(enhancers),
			ranges = IRanges::IRanges(start = enhancer_mids, end = enhancer_mids),
			name = GenomicRanges::mcols(enhancers)$name
		)
		dist_to_enh = GenomicRanges::distanceToNearest(mids_gr, enhancer_mids_gr)
		assigned_peaks$dist_to_enh = dist_to_enh@elementMetadata$distance
	}
	
	######################################################
	# Enrichment
	results = list()
	for (gobj in geneset_list) {
		message(sprintf("Test: Proximity test"))
		message(sprintf("Genesets: %s",gobj@type))
		message("Running tests..")
		
		if ("tss" %in% reglocation) {
			message("Running proximity to TSS test...")
			rtemp = test_proxReg(gobj, assigned_peaks, regloc = "tss", n_cores)
		} else if ("enhancer" %in% reglocation) {
			message("Running proximity to enhancers test...")
			rtemp = test_proxReg(gobj, assigned_peaks, regloc = "enhancer", n_cores)
		}
		
		# Annotate with geneset descriptions.
		rtemp$"Description" = as.character(mget(rtemp$Geneset.ID, gobj@set.name, ifnotfound=NA))
		rtemp$"Geneset.Type" = gobj@type
		
		results[[gobj@type]] = rtemp
	}
	enrich = Reduce(rbind,results)
	
	######################################################
	# Post-process enrichment
	# Order columns, add enriched/depleted column as needed, remove bad tests,
	# sort by p-value, rename rownames to integers
	enrich = post_process_enrichments(enrich)
	
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
		
		if (qc_plots) {
			filename_qcplots = file.path(out_path, sprintf("%s_qcplots.png", out_name))
			grDevices::png(filename_qcplots)
			#print(..plot_proxReg_spline(peaks = assigned_peaks, num_peaks = num_peaks))
			print(..plot_dist_to_tss(peakobj, tss))
			grDevices::dev.off()
			message("Wrote QC plots to: ",filename_qcplots)
		}
	}
	
	######################################################
	# Return objects as list
	return(list(
		peaks = assigned_peaks,
		results = enrich,
		opts = opts
	))
	
}