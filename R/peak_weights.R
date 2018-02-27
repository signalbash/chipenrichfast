
# This function calculates a weight for each peak according to how many loci it intersects.
calc_peak_weights = function(assigned_peaks, weighting) {
    assigned_peaks$peak_weights = 1
    if ("multiAssign" %in% weighting) {
        peak_counts = table(assigned_peaks$peak_id)
        peak_weights = data.frame(
            'peak_id' = names(peak_counts),
            'peak_weights' = as.numeric(1 / peak_counts),
            stringsAsFactors=FALSE)
	
        assigned_peaks = merge(assigned_peaks, peak_weights, by='peak_id')
        message("Assigning weights based on multiple peak assignments...")

    }
    
    if (any(c("signalValue", "logsignalValue") %in% weighting)) {
        if (!all(is.numeric(assigned_peaks$signalValue))) {
            stop("signalValue column is not all numeric!")
        } else if (any(assigned_peaks$signalValue<=0)) {
            stop("signalValue column should be all positive!")
        }
    }
    
    
    if ("signalValue" %in% weighting) {
        #assigned_peaks$peak_weights = assigned_peaks$peak_weights*log(assigned_peaks$signalValue)
        assigned_peaks$peak_weights = assigned_peaks$peak_weights*assigned_peaks$signalValue
        message("Assigning weights based on signalValue...")

    }
    if ("logsignalValue" %in% weighting) {
        #assigned_peaks$peak_weights = assigned_peaks$peak_weights*log(assigned_peaks$signalValue)
        assigned_peaks$peak_weights = assigned_peaks$peak_weights*log(assigned_peaks$signalValue)
        message("Assigning weights based on logsignalValue...")
        
    }
	return(assigned_peaks)
}


# Adds gene_weight column to the peaks per gene result
# Currently the gene_weight is the sum of the peak weights overlapping the locus
calc_genes_peak_weight = function(assigned_peaks, ppg) {
    message("Summing the peak weights...")
	# Sum up the normalized peak weights for each peak in a gene
    assigned_peaks$peak_weights = assigned_peaks$peak_weights/mean(assigned_peaks$peak_weights, na.rm=TRUE)
	rpg = stats::aggregate(peak_weights ~ gene_id, assigned_peaks, sum)
	
	d_rpg = data.frame(gene_id = rpg$gene_id, sum_peak_weight = rpg$peak_weight, stringsAsFactors=FALSE)
	
	result = merge(ppg, d_rpg, by='gene_id', all.x=TRUE)
	result$sum_peak_weight[is.na(result$sum_peak_weight)] = 0
	
	# Order by number of peaks in a gene.
	result = result[order(result$num_peaks,decreasing=TRUE),]
	
	return(result)
}