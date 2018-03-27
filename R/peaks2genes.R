# This function is used to create the *_peaks and *_peaks-per-gene files
# in case only the enrichment step has a problem. This way one does not
# need to rerun the making of these files whenever one just wants to test
# a new enrichment method.


#' @export
#' @include constants.R utils.R supported.R setup.R randomize.R
#' @include read.R assign_peaks.R peaks_per_gene.R
#' @include plot_dist_to_tss.R plot_chipenrich_spline.R

peaks2genes <- function(peaks,
                        out_name = "readyToEnrich",
                        out_path = getwd(),
                        genome = supported_genomes(),
                        genesets = c(
                            'GOBP',
                            'GOCC',
                            'GOMF'),
                        locusdef = "nearest_tss",
                        weighting = NULL,
                        mappability = NULL,
                        qc_plots = TRUE,
                        min_geneset_size = 15,
                        max_geneset_size = 2000,
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
