
# This function calculates a weight for each peak according to how many loci it intersects.
calc_peak_weights = function(assigned_peaks, weighting) {
    if (multiplicity %in% weighting) {
        peak_counts = table(assigned_peaks$peak_id)
        peak_weights = data.frame(
            'peak_id' = names(peak_counts),
            'peak_weight' = as.numeric(1 / peak_counts),
            stringsAsFactors=F)
	
        assigned_peaks = merge(assigned_peaks, peak_weights, by='peak_id')
    }
    
    if ("signalValue" %in% weighting) {
        assigned_peaks$peak_weights = assigned_peaks$peak_weights*log(assigned_peaks$signalValue)
    }
	return(assigned_peaks)
}


# Adds gene_weight column to the peaks per gene result
# Currently the gene_weight is the sum of the peak weights overlapping the locus
calc_genes_peak_weight = function(assigned_peaks, ppg) {
    message("Summing the peak weights...")
	# Sum up the peak weights for each peak in a gene
    logsignalValue = log(assigned_peaks$signalValue)
    assigned_peaks$peak_weight = logsignalValue/mean(logsignalValue, na.rm=T)
	rpg = stats::aggregate(peak_weight ~ gene_id, assigned_peaks, sum)
	
	d_rpg = data.frame(gene_id = rpg$gene_id, sum_peak_weight = rpg$peak_weight, stringsAsFactors=F)
	
	result = merge(ppg, d_rpg, by='gene_id', all.x=T)
	result$sum_peak_weight[is.na(result$sum_peak_weight)] = 0
	
	# Order by number of peaks in a gene.
	result = result[order(result$num_peaks,decreasing=T),]
	
	return(result)
}