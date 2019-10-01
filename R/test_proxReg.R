test_proxReg = function(geneset, peaks, regloc, n_cores) {
	# Restrict our genes/weights/peaks to only those genes in the genesets.
	# Here, geneset is not all combined, but GOBP, GOCC, etc.
	# i.e. A specific one.
	peaks = subset(peaks, peaks$gene_id %in% geneset@all.genes)
	
	# Run tests. NOTE: If os == 'Windows', n_cores is reset to 1 for this to work
	results_list = parallel::mclapply(as.list(ls(geneset@set.gene)), function(go_id) {
		single_proxReg(go_id, geneset, peaks, regloc)}, mc.cores = n_cores)
	
	# Collapse results into one table
	results = Reduce(rbind,results_list)
	
	# Correct for multiple testing
	results$FDR = p.adjust(results$P.value, method="BH")
	
	# Create enriched/depleted status column
	results$Status = ifelse(results$Effect > 0, 'closer', 'farther')
	
	results = results[order(results$P.value),]
	
	return(results)
}

single_proxReg = function(go_id, geneset, peaks, regloc) {

	# Genes in the geneset
	go_genes = geneset@set.gene[[go_id]]
	
	# Background genes and the background presence of a peak
	b_genes = peaks$gene_id %in% go_genes
	if (regloc == "tss"){
		x = peaks$scaled_dtss[b_genes]
		y = peaks$scaled_dtss[!b_genes]
	}
	if (regloc == "enhancer"){
		x = peaks$scaled_denh[b_genes]
		y = peaks$scaled_denh[!b_genes]
	}
	x <- x[is.finite(x)]
	y <- y[is.finite(y)]
	
	# Information about the geneset
	r_go_id = go_id
	r_go_genes_num = length(go_genes)

	# Information about peak genes
	go_genes_peak = peaks$gene_id[b_genes]
	r_go_genes_peak = length(table(go_genes_peak))
	r_go_genes_peak_num = length(go_genes_peak)

	r_effect = NA
	r_pval = NA
	
	if (length(x)== 0L | length(y) == 0L) {
		sprintf("Geneset: %s has zero peaks in one group. NAs given", go_id)
	} else {
		tryCatch(
			{ #Code from wilcox.test.default
			r <- rank(c(x, y))
			n.x <- as.double(length(x))
			n.y <- as.double(length(y))
			STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x + 1)/2)
			NTIES <- table(r)
			z <- STATISTIC - n.x * n.y/2
			SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) - 
							sum(NTIES^3 - NTIES)/((n.x + n.y) * (n.x + n.y - 
								1))))
			z = z/SIGMA
			
			# Results from the Wilcoxon test
			r_effect = -z; #We want closer to be positive effect
            r_pval = 2 * min(stats::pnorm(z),
                             stats::pnorm(z, lower.tail = FALSE))
			},
			error = {function(e) {warning(
				sprintf("Error in geneset: %s. NAs given", go_id))
			}}
		)
	}
	
	out = data.frame(
		"P.value"=r_pval,
		"Geneset ID"=r_go_id,
		"N Geneset Genes"=r_go_genes_num,
		"Geneset Peak Genes"=r_go_genes_peak,
		"N Geneset Peak Genes"=r_go_genes_peak_num,
		"Effect"=r_effect,
		stringsAsFactors = FALSE)
	
	return(out)
}
