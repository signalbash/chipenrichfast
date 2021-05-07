#' Run Proximity Regulation test on a set of narrow genomic regions
#' 
#' This method is designed for a set of narrow genomic regions (e.g. TF peaks) and is used to test
#' whether the genomic regions assigned to genes in a gene set are closer to 
#' regulatory locations (i.e. promoters or enhancers) than by chance.
#' 
#' @section Regulatory locations:
#' Current supported regulatory locations are gene transcription
#' start sites (tss) or enhancer locations (hg19 only)
#' 
#' @section Method:
#' ProxReg first calculates the distance between each peak midpoint and 
#' regulatory location in base pairs. For gene transcription start sites, 
#' since parts of the chromosome are more sparse than others, there is an
#' association with gene locus length that needs to be adjusted for.
#' When using tss as the regulatory location, the peak distances are 
#' adjusted for this confounding variable based on an average of 90 ENCODE
#' ChIP-seq experiments (details in citation pending). Similarly, for enhancers, 
#' distances depend on the density of enhancers within a gene locus, so distance
#' to enhancer is adjusted using an empirical average of 90 ChIP-seq ENCODE
#' experiments.
#' 
#' For each gene set of interest, the genomic regions are divided into two groups indicating
#' the gene with the nearest tss is in the gene set or not. A Wilcoxon Rank-Sum test is 
#' then done to test for a difference in the adjusted distances (either to tss or enhancer).
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
#' @param out_name Prefix string to use for naming output files. This should not
#' contain any characters that would be illegal for the system being used (Unix,
#' Windows, etc.) The default value is "proxReg", and a file "proxReg_results.tab"
#' is produced. If \code{qc_plots} is set, then a file "proxReg_qcplots.pdf"
#' is produced containing a number of quality control plots. If \code{out_name}
#' is set to NULL, no files are written, and results then must be retrieved from
#' the list returned by \code{proxReg}.
#' @param out_path Directory to which results files will be written out. Defaults
#' to the current working directory as returned by \code{\link{getwd}}.
#' @param genome One of the \code{supported_genomes()}. If reglocation = enhancer,
#' genome MUST be 'hg19'.
#' @param genesets A character vector of geneset databases to be tested for
#' enrichment. See \code{supported_genesets()}. Alternately, a file path to
#' a tab-delimited text file with header and first column being the geneset ID
#' or name, and the second column being Entrez Gene IDs. For an example custom
#' gene set file, see the vignette.
#' @param randomization One of: 'shuffle', 'unif', 'bylength', 'byenh'. These were used to
#' test for Type I error under the null hypothesis. A general user will never have to
#' use these.
#' @param reglocation One of: 'tss', 'enhancer'. Details in the "Regulatory locations" section
#' @param qc_plots A logical variable that enables the automatic generation of
#' plots for quality control.
#' @param min_geneset_size Sets the minimum number of genes a gene set may have
#' to be considered for testing.
#' @param max_geneset_size Sets the maximum number of genes a gene set may have
#' to be considered for testing.
#' @param n_cores The number of cores to use for testing. We recommend
#' using only up to the maximum number of \emph{physical} cores present, as
#' virtual cores do not significantly decrease runtime. Default number of cores
#' is set to 1. NOTE: Windows does not support multicore testing.
#' 
#' @return A list, containing the following items:
#'
#' \item{opts }{A data frame containing the arguments/values passed to \code{polyenrich}.}
#'
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
#'   \item{peak_id}{an ID given to unique combinations of chromosome, peak start, and peak end. }
#'   \item{chr}{the chromosome the peak originated from. }
#'   \item{peak_start}{start position of the peak. }
#'   \item{peak_end}{end position of the peak. }
#'   \item{gene_id}{the Entrez ID of the gene to which the peak was assigned. }
#'   \item{gene_symbol}{the official gene symbol for the gene_id (above). }
#'   \item{gene_locus_start}{the start position of the locus for the gene to which the peak was assigned (specified by the locus definition used.) }
#'   \item{gene_locus_end}{the end position of the locus for the gene to which the peak was assigned (specified by the locus definition used.) }
#'   \item{nearest_tss}{the closest TSS to this peak (for any gene, not necessarily the gene this peak was assigned to.) }
#'   \item{dist_to_tss}{the distance in bp to the closest TSS to this peak. }
#'   \item{nearest_tss_gene}{the gene having the closest TSS to the peak (should be the same as gene_id when using the nearest TSS locus definition.) }
#'   \item{nearest_tss_gene_strand}{the strand of the gene with the closest TSS. }
#'   \item{log_dtss}{ log of dist_to_tss}
#'   \item{log_gene_ll}{the log of length of the gene locus in bp }
#'   \item{scaled_dtss}{the adjusted distance to TSS, used in the calculations. Shown if reglocation = "tss" }
#'   \item{dist_to_enh}{the distance to the nearest enhancer. Shown if reglocation = "enhancer" }
#'   \item{avg_denh}{the empirical average for distance to the nearest enhancer for the gene the peak is assigned to. Shown if reglocation = enhancer }
#'   \item{scaled_denh}{the adjusted distance to the nearest enhancer. Shown if reglocation = enhancer }

#' }}
#'
#' \item{results }{
#' A data frame of the results from performing the proxReg test on
#' each geneset that was requested (all genesets are merged into one final data
#' frame.) The columns are:
#'
#' \describe{
#'   \item{Geneset.ID}{the identifier for a given gene set from the selected database.  For example, GO:0000003. }
#'   \item{Geneset.Type}{ specifies from which database the Geneset.ID originates.  For example, "Gene Ontology Biological Process."}
#'   \item{Description}{ gives a definition of the geneset. For example, "reproduction."}
#'   \item{P.Value}{the probability of observing the proxmity of genomic regions in the gene set given the null hypothesis that peaks are not closer or farther in the gene set.}
#'   \item{FDR}{the false discovery rate proposed by Bejamini \& Hochberg for adjusting the p-value to control for family-wise error rate.}
#'   \item{Effect}{the signed Wilcoxon statistic, with positive values meaning the gene set has closer genomic regions than expected by chance.}
#'   \item{Status}{specifies if the peaks in the gene set tend to be closer or farther than those not in the gene set.}
#'   \item{Odds.Ratio}{the estimated odds that peaks are associated with a given gene set compared to the odds that peaks are associated with other gene sets, after controlling for locus length and/or mappability.  An odds ratio greater than 1 indicates enrichment, and less than 1 indicates depletion.}
#'   \item{N.Geneset.Genes}{the number of genes in the gene set.}
#'   \item{N.Geneset.Peak.Genes}{the number of genes in the genes set that were assigned at least one peak.}
#'   \item{Geneset.Peak.Genes}{the list of genes from the gene set that had at least one peak assigned.}
#' }}
#' 
#' @examples 
#' # Run proxReg using an example dataset, assigning peaks to the nearest TSS,
#' # and on a small custom geneset
#' data(peaks_E2F4, package = 'chipenrich.data')
#' peaks_E2F4 = subset(peaks_E2F4, peaks_E2F4$chrom == 'chr1')
#' gs_path = system.file('extdata','vignette_genesets.txt', package='chipenrich')
#' results = proxReg(peaks_E2F4, reglocation = 'tss',
#' 			genome = 'hg19', genesets=gs_path, out_name=NULL)
#'
#' # Get the list of peaks that were assigned to genes and their distances to 
#' # regulatory regions.
#' assigned_peaks = results$peaks
#'
#' # Get the results of enrichment testing.
#' enrich = results$results
#' 
#' @export
#' @include constants.R utils.R supported.R setup.R
#' @include read.R assign_peaks.R 
#' @include plot_dist_to_tss.R
#' @include test_proxReg.R
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
	randomization = NULL,
	qc_plots = TRUE,
	min_geneset_size = 15,
	max_geneset_size = 2000,
	n_cores = 1
) {
	# Kill if reglocation is not tss or enhancer
	if (!(reglocation %in% c("tss","enhancer"))) {
		stop("Unsupported regulatory location!")
	}
	
	# Set genome to hg due to enhancer-specific 
	if ("enhancer" %in% reglocation & genome != "hg19" & genome != "hg38") {
		stop("Genome must be hg19 or hg38 to use enhancer regulatory location")
	}
	
	genome = match.arg(genome)
	
	n_cores = reset_ncores_for_windows(n_cores)
	
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
	
	#Stop if enhancer reglocation is selected and genome is not hg19
	if (reglocation == "enhancer" & genome != "hg19" & genome != "hg38") {
		stop("Enhancer regulatory location only compatible with hg19 or hg38 genome")
	}
	
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
		assigned_peaks$log_gene_ll = log(assigned_peaks$gene_locus_end-assigned_peaks$gene_locus_start)
		pred_log_dtss = as.numeric(mgcv::predict.gam(chipenrich.data::spline.log_dtss.90ENCODE, assigned_peaks, type="link"))
		assigned_peaks$scaled_dtss = assigned_peaks$log_dtss-pred_log_dtss
	}
	
	# Creating a column for enhancer distances
	if ("enhancer" %in% reglocation) {
		enhancers = chipenrich.data::enhancer.hg19
		gene.enh.desc = chipenrich.data::gene.enh.desc.hg19
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
		assigned_peaks$log_dtss = log(abs(assigned_peaks$dist_to_tss)+1)
		assigned_peaks$log_gene_ll = log(assigned_peaks$gene_locus_end-assigned_peaks$gene_locus_start)
		pred_log_dtss = as.numeric(mgcv::predict.gam(chipenrich.data::spline.log_dtss.90ENCODE, assigned_peaks, type="link"))
		assigned_peaks$avgdenh = sapply(assigned_peaks$gene_id, 
										function(x) {a = gene.enh.desc$avg_denh_emp[gene.enh.desc$gene_id == x]
                                            return(ifelse(length(a)==0,NA,a))
                                        })
        if (any(is.na(assigned_peaks$avgdenh))) { #Impute the missing genes as average of the others
        	dlm = lm(avgdenh ~ log_gene_ll, data = assigned_peaks)
            assigned_peaks$avgdenh[is.na(assigned_peaks$avgdenh)] = stats::predict.lm(dlm, assigned_peaks[is.na(assigned_peaks$avgdenh),])
     		rm(dlm)
        }
		assigned_peaks$scaled_denh = log(abs(assigned_peaks$dist_to_enh)+1) - log(assigned_peaks$avgdenh+1)
	}
	
	#Randomizations
	if (!is.null(randomization)){
		if (randomization == "shuffle") { #Just shuffle gene ids
			message("Randomizing by shuffling gene ids...")
			assigned_peaks$gene_id = sample(assigned_peaks$gene_id)
		} else if (randomization == "unif") { #Uniformly pick from random available genes
			message("Randomizing by uniformly picking gene ids...")
			assigned_peaks$gene_id = as.integer(sample(names(table(assigned_peaks$gene_id)), nrow(assigned_peaks), replace = T))
		} else if (randomization == "bylength") { #bins of locus length, style of unif
			message("Randomizing by picking gene ids by locus length...")
			peakstemp = assigned_peaks[,c("peak_id","gene_id","gene_locus_start","gene_locus_end")]
			peakstemp$gene_locus_length = peakstemp$gene_locus_end-peakstemp$gene_locus_start
			peakstemp$dupe_genes = duplicated(peakstemp$gene_id)
			peakstemp = peakstemp[sample(1:nrow(peakstemp)),]
			peakstemp2 = peakstemp[!peakstemp$dupe_genes,]
			peakstemp2 = peakstemp2[order(peakstemp2$gene_locus_length),]
			rownames(peakstemp2) = 1:nrow(peakstemp2)
			peakstemp2$group = floor((as.numeric(rownames(peakstemp2))+99)/100)
			peakstemp$group = sapply(peakstemp$gene_id, function(x){peakstemp2$group[peakstemp2$gene_id==x]})
			peakstemp$gene_id = sapply(peakstemp$group, function(x){sample(peakstemp2$gene_id[peakstemp2$group==x],1)})
			assigned_peaks$gene_id_pre = assigned_peaks$gene_id
			assigned_peaks = merge(assigned_peaks[,-5], peakstemp[,c("peak_id","gene_id")], by = "peak_id")
		} else if (randomization == "byenh"){
			message("Randomizing by picking gene ids by enhancer dist...")
			peakstemp = assigned_peaks[,c("peak_id","gene_id","avgdenh")]
			peakstemp$dupe_genes = duplicated(peakstemp$gene_id)
			peakstemp = peakstemp[sample(1:nrow(peakstemp)),]
			peakstemp2 = peakstemp[!peakstemp$dupe_genes,]
			peakstemp2 = peakstemp2[sample(1:nrow(peakstemp2)),]
			peakstemp2 = peakstemp2[order(peakstemp2$avgdenh),]
			rownames(peakstemp2) = 1:nrow(peakstemp2)
			peakstemp2$group = floor((as.numeric(rownames(peakstemp2))+499)/500)
			peakstemp$group = sapply(peakstemp$gene_id, function(x){peakstemp2$group[peakstemp2$gene_id==x]})
			peakstemp$gene_id = sapply(peakstemp$group, function(x){sample(peakstemp2$gene_id[peakstemp2$group==x],1)})
			assigned_peaks$gene_id_pre = assigned_peaks$gene_id
			assigned_peaks = merge(assigned_peaks[,-5], peakstemp[,c("peak_id","gene_id")], by = "peak_id")
			
		} else {
			stop("Unsupported randomization!")
		}
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
