#' Run the test process up to, but not including the enrichment tests.
#'
#' This function is used to create the *_peaks and *_peaks-per-gene files
#' This way one does not need to remake these files whenever one just wants
#' to test enrichment methods.
#'
#' @section Randomizations:
#' Randomization of locus definitions allows for the assessment of Type I Error
#' under the null hypothesis. The randomization codes are:
#' \describe{
#'	\item{\code{NULL}:}{ No randomizations, the default.}
#' 	\item{'complete':}{ Shuffle the \code{gene_id} and \code{symbol} columns of the
#' \code{locusdef} together, without regard for the chromosome location, or locus length.
#' The null hypothesis is that there is no true gene set enrichment.}
#' 	\item{'bylength':}{ Shuffle the \code{gene_id} and \code{symbol} columns of the
#' \code{locusdef} together within bins of 100 genes sorted by locus length. The null
#' hypothesis is that there is no true gene set enrichment, but with preserved locus
#' length relationship.}
#' 	\item{'bylocation':}{ Shuffle the \code{gene_id} and \code{symbol} columns of the
#' \code{locusdef} together within bins of 50 genes sorted by genomic location. The null
#' hypothesis is that there is no true gene set enrichment, but with preserved
#' genomic location.}
#' }
#' The return value with a selected randomization is the same list as without.
#' To assess the Type I error, the \code{alpha} level for the particular data set
#' can be calculated by dividing the total number of gene sets with p-value < \code{alpha}
#' by the total number of tests. Users may want to perform multiple randomizations
#' for a set of peaks and take the median of the \code{alpha} values.
#'
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
#' Windows, etc.) The default value is "polyenrich", and a file "polyenrich_results.tab"
#' is produced. If \code{qc_plots} is set, then a file "polyenrich_qcplots.pdf"
#' is produced containing a number of quality control plots. If \code{out_name}
#' is set to NULL, no files are written, and results then must be retrieved from
#' the list returned by \code{polyenrich}.
#' @param out_path Directory to which results files will be written out. Defaults
#' to the current working directory as returned by \code{\link{getwd}}.
#' @param genome One of the \code{supported_genomes()}.
#' @param locusdef One of: 'nearest_tss', 'nearest_gene', 'exon', 'intron', '1kb',
#' '1kb_outside', '1kb_outside_upstream', '5kb', '5kb_outside', '5kb_outside_upstream',
#' '10kb', '10kb_outside', '10kb_outside_upstream'. For a description of each,
#' see the vignette or \code{\link{supported_locusdefs}}. Alternately, a file path for
#' a custom locus definition. NOTE: Must be for a \code{supported_genome()}, and
#' must have columns 'chr', 'start', 'end', and 'gene_id' or 'geneid'. For an
#' example custom locus definition file, see the vignette.
#' @param weighting (Poly-Enrich only) character string specifying the weighting method if method is
#' chosen to be 'polyenrich_weighted'. Current options are: 'signalValue',
#' 'logsignalValue', and 'multiAssign'.
#' @param mappability One of \code{NULL}, a file path to a custom mappability file,
#' or an \code{integer} for a valid read length given by \code{supported_read_lengths}.
#' If a file, it should contain a header with two column named 'gene_id' and 'mappa'.
#' Gene IDs should be Entrez IDs, and mappability values should range from 0 and 1.
#' For an example custom mappability file, see the vignette. Default value is NULL.
#' @param qc_plots A logical variable that enables the automatic generation of
#' plots for quality control.
#' @param num_peak_threshold (ChIP-Enrich only) Sets the threshold for how many peaks a gene must
#' have to be considered as having a peak. Defaults to 1. Only relevant for
#' Fisher's exact test and ChIP-Enrich methods.
#' @param randomization One of \code{NULL}, 'complete', 'bylength', or 'bylocation'.
#' See the Randomizations section below.
#'
#' @section Poly-Enrich Weighting Options:
#' Poly-Enrich also allows weighting of individual peaks. Currently the options are:
#' \describe{
#'  \item{'signalValue:'}{ weighs each peak based on the Signal Value given in the
#' narrowPeak format or a user-supplied column, normalized to have mean 1.}
#'  \item{'logsignalValue:'}{ weighs each peak based on the log Signal Value given in the
#' narrowPeak format or a user-supplied column, normalized to have mean 1.}
#'  \item{'multiAssign:'}{ weighs each peak by the inverse of the number of genes
#' it is assigned to.}
#' }
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
#'   \item{peak_id}{ is an ID given to unique combinations of chromosome, peak start, and peak end. }
#'   \item{chr}{ is the chromosome the peak originated from. }
#'   \item{peak_start}{ is start position of the peak. }
#'   \item{peak_end}{ is end position of the peak. }
#'   \item{peak_midpoint}{ is the midpoint of the peak. }
#'   \item{gene_id}{ is the Entrez ID of the gene to which the peak was assigned. }
#'   \item{gene_symbol}{ is the official gene symbol for the gene_id (above). }
#'   \item{gene_locus_start}{ is the start position of the locus for the gene to which the peak was assigned (specified by the locus definition used.) }
#'   \item{gene_locus_end}{ is the end position of the locus for the gene to which the peak was assigned (specified by the locus definition used.) }
#'   \item{nearest_tss}{ is the closest TSS to this peak (for any gene, not necessarily the gene this peak was assigned to.) }
#'   \item{nearest_tss_gene}{ is the gene having the closest TSS to the peak (should be the same as gene_id when using the nearest TSS locus definition.) }
#'   \item{nearest_tss_gene_strand}{ is the strand of the gene with the closest TSS. }
#' }}
#'
#' \item{peaks_per_gene }{
#' A data frame of the count of peaks per gene. The columns are:
#'
#' \describe{
#'   \item{gene_id}{ is the Entrez Gene ID. }
#'   \item{length}{ is the length of the gene's locus (depending on which locus definition you chose.)}
#'   \item{log10_length}{ is the log10(locus length) for the gene.}
#'   \item{num_peaks}{ is the number of peaks that were assigned to the gene, given the current locus definition. }
#'   \item{peak}{ is whether or not the gene has a peak. }
#' }}
#'
#' @examples
#'
#' # Run peaks2genes using an example dataset, assigning peaks to the nearest TSS
#' data(peaks_E2F4, package = 'chipenrich.data')
#' peaks_E2F4 = subset(peaks_E2F4, peaks_E2F4$chrom == 'chr1')
#' gs_path = system.file('extdata', package='chipenrich')
#' results = peaks2genes(peaks_E2F4, locusdef='nearest_tss',
#' 			genome = 'hg19', out_name=NULL)
#'
#' # Get the list of peaks that were assigned to genes.
#' assigned_peaks = results$peaks
#'
#' @export
#' @include constants.R utils.R supported.R setup.R randomize.R
#' @include read.R assign_peaks.R peaks_per_gene.R
#' @include plot_dist_to_tss.R plot_chipenrich_spline.R

peaks2genes <- function(peaks,
                        out_name = "readyToEnrich",
                        out_path = getwd(),
                        genome = supported_genomes(),
                        locusdef = "nearest_tss",
                        weighting = NULL,
                        mappability = NULL,
                        qc_plots = TRUE,
                        num_peak_threshold = 1,
                        randomization = NULL
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
    # Setup locus definitions, and mappability

    ldef_list = setup_locusdef(locusdef, genome, randomization)
    ldef = ldef_list[['ldef']]
    tss = ldef_list[['tss']]

    mappa = setup_mappa(mappa_code = mappability, genome = genome, ldef_code = locusdef, ldef_obj = ldef)


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

    ######################################################
    # Compute peaks per gene table
    ppg = num_peaks_per_gene(assigned_peaks, ldef, mappa)

    if (!is.null(weighting)) {
        if (!all(weighting %in% c("logsignalValue","signalValue","multiAssign"))) {
            # Unsupported weights
            stop(sprintf("Unsupported weights: %s",
            paste(weighting[which(!(weighting %in% c("logsignalValue","signalValue","multiAssign")))],collapse=", ")))
        }
        assigned_peaks = calc_peak_weights(assigned_peaks, weighting)
        ppg = calc_genes_peak_weight(assigned_peaks, ppg)

    }

    ######################################################
    # Write result objects to files
    if (!is.null(out_name)) {
        filename_peaks = file.path(out_path, sprintf("%s_peaks.tab", out_name))
        write.table(assigned_peaks, file = filename_peaks, row.names = FALSE, quote = FALSE, sep = "\t")
        message("Wrote peak-to-gene assignments to: ", filename_peaks)

        filename_opts = file.path(out_path, sprintf("%s_opts.tab", out_name))
        write.table(opts, file = filename_opts, row.names = FALSE, quote = FALSE, sep = "\t")
        message("Wrote run options/arguments to: ", filename_opts)

        filename_ppg = file.path(out_path, sprintf("%s_peaks-per-gene.tab", out_name))
        write.table(ppg, file = filename_ppg, row.names = FALSE, quote = FALSE, sep = "\t")
        message("Wrote count of peaks per gene to: ", filename_ppg)

        if (qc_plots) {
            #filename_qcplots = file.path(out_path, sprintf("%s_qcplots.png", out_name))
            #grDevices::png(filename_qcplots)
            #print(..plot_chipenrich_spline(gpw = ppg, mappability = mappability, num_peaks = num_peaks))
            #print(..plot_dist_to_tss(peakobj, tss))
            #grDevices::dev.off()


            filename_qcplots_chip = file.path(out_path, sprintf("%s_qcplots_chip.png", out_name))
            filename_qcplots_poly = file.path(out_path, sprintf("%s_qcplots_poly.png", out_name))
            filename_disttotss = file.path(out_path,sprintf("%s_locuslength.png",out_name));


            grDevices::png(filename_qcplots_chip)
            print(..plot_chipenrich_spline(gpw = ppg, mappability = mappability, num_peaks = num_peaks))
            grDevices::dev.off()

            grDevices::png(filename_qcplots_poly)
            print(..plot_polyenrich_spline(gpw = ppg, mappability = mappability, num_peaks = num_peaks))
            grDevices::dev.off()

            grDevices::png(filename_disttotss)
            print(..plot_dist_to_tss(peakobj, tss))
            grDevices::dev.off()


            message("Wrote QC plots to: ",filename_qcplots_poly)
        }
    }

    ######################################################
    # Return objects as list
    return(list(
    peaks = assigned_peaks,
    opts = opts,
    peaks_per_gene = ppg
    ))
}
