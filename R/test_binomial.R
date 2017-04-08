test_binomial = function(geneset, ppg) {
	warning("The binomial test is provided for comparison purposes only. It will almost always give biased results favoring gene sets with short average locus length.")

	# Restrict to genes only within the genesets.
	ppg = subset(ppg, ppg$gene_id %in% geneset@all.genes)

	# Total # of peaks.
	total_peaks = sum(ppg$num_peaks)

	# Gene IDs must be character to be stored into an environment as names.
	ppg$gene_id = as.character(ppg$gene_id)

	# Index gene --> peak count.
	gene2peakc = new.env(parent = emptyenv())
	apply(ppg, 1, function(x) assign(x["gene_id"], x["num_peaks"], envir = gene2peakc))

	# Index gene --> log10 length.
	gene2length = new.env(parent = emptyenv())
	apply(ppg, 1, function(x) assign(x["gene_id"], x["length"], envir = gene2length))

	# Length of the genome.
	genome_length = sum(as.numeric(ppg$length))

	estimates = c()
	pvals = c()
	goids = c()
	gonums = c()
	n_go_peak_genes = c()
	r_go_peak_genes = c()
	for (i in 1:length(geneset@set.gene)) {
		go_term = ls(geneset@set.gene)[i]
		go_genes = as.character(geneset@set.gene[[go_term]])

		# Eliminate GO genes that aren't in the ppg.
		#go_genes = Filter(function(x) x %in% ppg$gene_id,go_genes)
		go_genes = go_genes[go_genes %in% ppg$gene_id]

		go_gene_lengths = as.numeric(sapply(go_genes,function(x) get(x, gene2length)))
		go_term_length = sum(go_gene_lengths)

		go_gene_counts = as.numeric(sapply(go_genes,function(x) get(x, gene2peakc)))
		go_term_count = sum(go_gene_counts)

		binom_p = go_term_length / genome_length
		btest = stats::binom.test(go_term_count, total_peaks, binom_p, alternative = "greater")

		goids[i] = go_term
		pvals[i] = btest$p.value
		estimates[i] = btest$estimate
		gonums[i] = length(go_genes)

		go_peak_genes = ppg[(ppg$gene_id %in% go_genes) & (ppg$num_peaks > 0),]$gene_id
		r_go_peak_genes[i] = paste(go_peak_genes, collapse = ", ")
		n_go_peak_genes[i] = length(go_peak_genes)
	}

	fdr = stats::p.adjust(pvals,method="BH")

	results = data.frame(
		"Geneset ID" = goids,
		"N Geneset Genes" = gonums,
		"N Geneset Peak Genes" = n_go_peak_genes,
		"P.Success" = estimates,
		"P-value" = pvals,
		"FDR" = fdr,
		"Geneset Peak Genes" = r_go_peak_genes,
		stringsAsFactors = FALSE)

	results = results[order(results$P.value), ]
	return(results)
}
