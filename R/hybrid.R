#' Running Hybrid test, either from scratch or using two results files
#'
#' Hybrid test is designed for people unsure of which test between ChIP-Enrich
#' and Poly-Enrich to use, so it takes information of both and gives adjusted
#' P-values. For more about ChIP- and Poly-Enrich, consult their corresponding
#' documentation.
#'
#' @section Hybrid p-values:
#' Given n tests that test for the same hypothesis, same Type I error rate, and
#' converted to p-values: \code{p_1, ..., p_n}, the Hybrid p-value is computed as:
#' \code{n*min(p_1, ..., p_n)}. This hybrid test will have at most the same
#' Type I error as any individual test, and if any of the tests have 100\% power as
#' sample size goes to infinity, then so will the hybrid test.
#'
#' @section Function inputs:
#' Every input in hybridenrich is the same as in chipenrich and polyenrich. Inputs
#' unique to chipenrich are: num_peak_threshold; and inputs unique to polyenrich are:
#' weighting. Currently the test only supports running chipenrich and polyenrich, but
#' future plans will allow you to run any number of different support tests.
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
#' or name, and the second column being Entrez Gene IDs. For an example custom
#' gene set file, see the vignette.
#' @param locusdef One of: 'nearest_tss', 'nearest_gene', 'exon', 'intron', '1kb',
#' '1kb_outside', '1kb_outside_upstream', '5kb', '5kb_outside', '5kb_outside_upstream',
#' '10kb', '10kb_outside', '10kb_outside_upstream'. For a description of each,
#' see the vignette or \code{\link{supported_locusdefs}}. Alternately, a file path for
#' a custom locus definition. NOTE: Must be for a \code{supported_genome()}, and
#' must have columns 'chr', 'start', 'end', and 'gene_id' or 'geneid'. For an
#' example custom locus definition file, see the vignette.
#' @param methods A character string array specifying the method to use for enrichment
#' testing. Currently actually unused as the methods are forced to be one chipenrich
#' and one polyenrich.
#' @param weighting A character string specifying the weighting method. Method name will
#' automatically be "polyenrich_weighted" if given weight options. Current options are:
#' 'signalValue', 'logsignalValue', and 'multiAssign'.
#' @param mappability One of \code{NULL}, a file path to a custom mappability file,
#' or an \code{integer} for a valid read length given by \code{supported_read_lengths}.
#' If a file, it should contain a header with two column named 'gene_id' and 'mappa'.
#' Gene IDs should be Entrez IDs, and mappability values should range from 0 and 1.
#' For an example custom mappability file, see the vignette. Default value is NULL.
#' @param qc_plots A logical variable that enables the automatic generation of
#' plots for quality control.
#' @param min_geneset_size Sets the minimum number of genes a gene set may have
#' to be considered for enrichment testing.
#' @param max_geneset_size Sets the maximum number of genes a gene set may have
#' to be considered for enrichment testing.
#' @param num_peak_threshold Sets the threshold for how many peaks a gene must
#' have to be considered as having a peak. Defaults to 1. Only relevant for
#' Fisher's exact test and ChIP-Enrich methods.
#' @param randomization One of \code{NULL}, 'complete', 'bylength', or 'bylocation'.
#' See the Randomizations section below.
#' @param n_cores The number of cores to use for enrichment testing. We recommend
#' using only up to the maximum number of \emph{physical} cores present, as
#' virtual cores do not significantly decrease runtime. Default number of cores
#' is set to 1. NOTE: Windows does not support multicore enrichment.
#'
#' @section Joining two results files:
#' Combines two existing results files and returns one results file with hybrid
#' p-values and FDR included. Current allowed inputs are objects from any of
#' the supplied enrichment tests or a dataframe with at least the following columns:
#' \code{P.value, Geneset.ID}. Optional columns include: \code{Status}. Currently
#' we only allow for joining two results files, but future plans will allow you to join
#' any number of results files.
#'
#' @return A data.frame containing:
#' \item{results }{
#' A data frame of the results from performing the gene set enrichment test on
#' each geneset that was requested (all genesets are merged into one final data
#' frame.) The columns are:
#'
#' \describe{
#'   \item{Geneset.ID}{ is the identifier for a given gene set from the selected database.  For example, GO:0000003. }
#'   \item{P.Value.x}{ is the probability of observing the degree of enrichment of the gene set given the null hypothesis
#'                     that peaks are not associated with any gene sets, for the first test.}
#'   \item{P.Value.y}{ is the same as above except for the second test.}
#'   \item{P.Value.Hybrid}{ The calculated Hybrid p-value from the two tests}
#'   \item{FDR.Hybrid}{ is the false discovery rate proposed by Bejamini \& Hochberg for adjusting the p-value to control for family-wise error rate.}
#'  Other variables given will also be included, see the corresponding methods' documentation for their details.
#' }}
#'
#' @export
#' @include chipenrich.R polyenrich.R

hybridenrich <- function(	peaks,
    out_name = "hybridenrich",
    out_path = getwd(),
    genome = supported_genomes(),
    genesets = c('GOBP','GOCC','GOMF'),
    locusdef = "nearest_tss",
    methods = c('chipenrich','polyenrich'),
    weighting = NULL,
    mappability = NULL,
    qc_plots = TRUE,
    min_geneset_size = 15,
    max_geneset_size = 2000,
    num_peak_threshold = 1,
    randomization = NULL,
    n_cores = 1
) {
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
    
    ldef_list = setup_locusdef(locusdef, genome, randomization)
    ldef = ldef_list[['ldef']]
    tss = ldef_list[['tss']]
    
    geneset_list = setup_genesets(gs_codes = genesets, ldef_obj = ldef, genome = genome, min_geneset_size = min_geneset_size, max_geneset_size = max_geneset_size)
    
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
    # There is a for loop here to run through more than two tests, but later code only supports CE and PE for now. You need
    # to change the hybrid.join lines before deciding to add hybrid functionality to more than 2 methods!
    methodindex = 0
    
    # Initializing to be a list later.
    enrich = NULL
    for (method in methods) {
        methodindex = methodindex+1
        if (!is.null(weighting) & method=="polyenrich") {
            warning("Weights given but polyenrich_weighted method not chosen. Automatically switching to weighted.")
            method = "polyenrich_weighted"
        } # This isn't exactly the best workaround, as it forces weighted even if you don't want to. We need to change this
        ### if we want to hybrid between weighted and non-weighted. We also cannot hybrid two different types of weighting.
   
        ############################################################################
        # CHECK method and get() it if okay
        testf = get_test_method(method)
        test_func = get(testf)
        method_name = METHOD_NAMES[[method]]

        ######################################################
        # Enrichment
        results = list()
        for (gobj in geneset_list) {
            message(sprintf("Test: %s",method_name))
            message(sprintf("Genesets: %s",gobj@type))
            message("Running tests..")
            if (testf == "test_chipenrich_slow") {
                rtemp = test_func(gobj,ppg,n_cores)
            }
            if (testf == "test_fisher_exact") {
                rtemp = test_func(gobj,ppg,alternative="two.sided")
            }
            if (testf == "test_binomial") {
                rtemp = test_func(gobj,ppg)
            }
            if (testf == "test_approx") {
                rtemp = test_func(gobj,ppg,nwp=FALSE,n_cores)
            }
            if (testf == "test_chipenrich") {
                rtemp = test_func(gobj,ppg,n_cores)
            }
            if (testf == "test_chipapprox") {
                rtemp = test_func(gobj,ppg,n_cores)
            }
            if (testf == "test_polyenrich_slow") {
                rtemp = test_func(gobj,ppg,n_cores)
            }
            if (testf == "test_polyenrich") {
                rtemp = test_func(gobj,ppg,n_cores)
            }
            if (testf == "test_polyenrich_weighted") {
                # Note that there is an extra input for the name of the count column
                # In the future when we make a function purely for enrichment, this will be where you
                # give the name of the count column.
                rtemp = test_func(gobj,ppg,n_cores, "sum_peak_weight")
            }
            if (testf == "test_polyapprox") {
                rtemp = test_func(gobj,ppg,n_cores)
            }

            
            # Annotate with geneset descriptions.
            rtemp$"Description" = as.character(mget(rtemp$Geneset.ID, gobj@set.name, ifnotfound=NA))
            rtemp$"Geneset.Type" = gobj@type
            
            results[[gobj@type]] = rtemp
        }
        enrich[[methodindex]] = Reduce(rbind,results)
        
        ######################################################
        # Post-process enrichment
        # Order columns, add enriched/depleted column as needed, remove bad tests,
        # sort by p-value, rename rownames to integers
        enrich[[methodindex]] = post_process_enrichments(enrich[[methodindex]])
    
    }
    
    #Temporary only 2 tests hybridenrich solution:
    message("Combining tests..")
    hybrid = hybrid.join(enrich[[1]],enrich[[2]])
    
    
    ######################################################
    # Write result objects to files
    if (!is.null(out_name)) {
        filename_analysis = file.path(out_path, sprintf("%s_results.tab", out_name))
        write.table(hybrid, file = filename_analysis, row.names = FALSE, quote = FALSE, sep = "\t")
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
        
        if (qc_plots) {
            filename_qcplots = file.path(out_path, sprintf("%s_qcplots.png", out_name))
            #grDevices::png(filename_qcplots)
            #print(..plot_chipenrich_spline(gpw = ppg, mappability = mappability, num_peaks = num_peaks))
            #print(..plot_dist_to_tss(peakobj, tss))
            #grDevices::dev.off()
            
            
            filename_qcplots_chip = file.path(out_path, sprintf("%s_qcplots_chip.png", out_name))
            filename_qcplots_poly = file.path(out_path, sprintf("%s_qcplots_poly.png", out_name))
            filename_disttotss = file.path(out_path,sprintf("%s_locuslength.jpeg",out_name));
            
            
            grDevices::png(filename_qcplots_chip)
            print(..plot_chipenrich_spline(gpw = ppg, mappability = mappability, num_peaks = num_peaks))				
            grDevices::dev.off()
            
            grDevices::png(filename_qcplots_poly)
            print(..plot_polyenrich_spline(gpw = ppg, mappability = mappability, num_peaks = num_peaks))
            grDevices::dev.off()
            
            grDevices::png(filename_disttotss)
            print(..plot_dist_to_tss(peakobj, tss))
            grDevices::dev.off()
            
            
            message("Wrote QC plots to: ",filename_qcplots)
        }
    }
    
    ######################################################
    # Return objects as list
    return(list(
    peaks = assigned_peaks,
    results = hybrid,
    opts = opts,
    peaks_per_gene = ppg
    ))
}


hybridenrich.old <- function(	peaks,
						out_name = "hybridenrich",
						out_path = getwd(),
						genome = supported_genomes(),
						genesets = c(
							'GOBP',
							'GOCC',
							'GOMF'),
						locusdef = "nearest_tss",
						methods = c('chipenrich','polyenrich'),
						weighting = NULL,
						mappability = NULL,
						qc_plots = TRUE,
						min_geneset_size = 15,
						max_geneset_size = 2000,
						num_peak_threshold = 1,
						randomization = NULL,
						n_cores = 1
) {
    if (!is.null(out_name)) {
        out_chip = sprintf("%s_chip",out_name)
        out_poly = sprintf("%s_poly",out_name)
    } else {
        out_chip = NULL
        out_poly = NULL
    }
    
    #Check if methods number isn't 2. Will support more than 2 later.
    #Also not really relevant as you're forced to use chip and poly anyway.
    if (length(methods) != 2) {
        stop("Hybrid test currently only supports exactly two methods!")
    }
    
	results1 = chipenrich(
        peaks = peaks,
        out_name = out_chip,
        out_path = out_path,
        genome = genome,
        genesets = genesets,
        locusdef = locusdef,
        method = "chipenrich",
        mappability = mappability,
        qc_plots = qc_plots,
        min_geneset_size = min_geneset_size,
        max_geneset_size = max_geneset_size,
        num_peak_threshold = num_peak_threshold,
        randomization = randomization,
        n_cores = n_cores)
        
        
    if (is.null(weighting)) {
        polymeth = "polyenrich"
    } else {
        polymeth = "polyenrich_weighted"
    }
    results2 = polyenrich(
        peaks = peaks,
        out_name = out_poly,
        out_path = out_path,
        genome = genome,
        genesets = genesets,
        locusdef = locusdef,
        method = polymeth,
        weighting = weighting,
        mappability = mappability,
        qc_plots = qc_plots,
        min_geneset_size = min_geneset_size,
        max_geneset_size = max_geneset_size,
        randomization = randomization,
        n_cores = n_cores)
        
    hybrid = hybrid.join(results1,results2)
    
    return(hybrid)
}



#User gives results objects (object or just results), and appends hybrid results
hybrid.join <- function(test1, test2) {
	#Check if they inputed the test object or just the results file, checked by seeing
    # if object has a $results part
    if ("results" %in% names(test1)) {
        #If entire object, extract the results section
        results1 = test1$results
    } else if ("P.value" %in% names(test1)) {
        results1 = test1
    } else {
        stop("First object is not a valid output or does not have P.value column")
    }
	
    if ("results" %in% names(test2)) {
        #If entire object, extract the results section
        results2 = test2$results
    } else if ("P.value" %in% names(test2)) {
        results2 = test2
    } else {
        stop("Second object is not a valid output or does not have P.value column")
    }
    
    #Check for Geneset.ID column
    if (!("Geneset.ID" %in% names(results1))) {
        stop("First object does not have Geneset.ID column")
    }
    if (!("Geneset.ID" %in% names(results2))) {
        stop("Second object does not have Geneset.ID column")
    }
    
    
    

    
    #Separate tree if the data does not have Status column
    if ("Status" %in% names(results1) & "Status" %in% names(results2)) {
        #Extract p-value and status of both tests
        Pvals1 = results1[,c("Geneset.ID","P.value","Status")]
        Pvals2 = results2[,c("Geneset.ID","P.value","Status")]
    
    } else {
        #Extract p-value only
        Pvals1 = results1[,c("Geneset.ID","P.value")]
        Pvals2 = results2[,c("Geneset.ID","P.value")]

    }

    #Merge by Geneset.ID
    PvalsH = merge(Pvals1, Pvals2, by="Geneset.ID")
    #If 0 remain, stop.
    if (nrow(PvalsH) == 0) {
        stop("No common genesets in the two datasets!")
    }
    message(sprintf("Total of %s common Geneset.IDs", nrow(PvalsH)))


    PvalsH$P.value.Hybrid = 2*pmin(PvalsH$P.value.x, PvalsH$P.value.y, 0.5)
	
	#Run B-H to adjust for FDR for hybrid p-values
    PvalsH$FDR.Hybrid = stats::p.adjust(PvalsH$P.value.Hybrid, method = "BH")
    
    
    #Include enrich/depleted status if available and combine
    if ("Status" %in% names(results1) & "Status" %in% names(results2)) {
        PvalsH$Status.Hybrid = ifelse(PvalsH$Status.x == PvalsH$Status.y, PvalsH$Status.x, "Inconsistent")
        #Combine both results files together and append hybrid p-value and FDR
        resultsH = merge(results1[,-which(colnames(results1) %in% c("P.value","Status"))], PvalsH, by = "Geneset.ID")
    } else {
        resultsH = merge(results1[,-which(colnames(results1) %in% c("P.value"))], PvalsH, by = "Geneset.ID")
    }
    
    
    #Reorder columns?????
	
	#Output final results
	return(resultsH)
	
}