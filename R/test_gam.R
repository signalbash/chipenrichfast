# Extracted contents of for loop implementation in test_gam and test_gam_ratio to
# enable multicore processing of genesets. This function gets called in mclapply
# in both test_gam and test_gam_ratio, while using the correct model, and reporting
# the correct information based on the method parameter.
single_gam = function(go_id, geneset, gpw, method, model) {
	final_model = as.formula(model)

	# Genes in the geneset
	go_genes = geneset@set.gene[[go_id]]

	# Background genes and the background presence of a peak
	b_genes = gpw$gene_id %in% go_genes
	sg_go = gpw$peak[b_genes]

 	# Information about the geneset
	r_go_id = go_id
	r_go_genes_num = length(go_genes)
	r_go_genes_avg_length = mean(gpw$length[b_genes])

    # Information about peak genes
	go_genes_peak = gpw$gene_id[b_genes][sg_go == 1]
	r_go_genes_peak = paste(go_genes_peak, collapse = ", ")
	r_go_genes_peak_num = length(go_genes_peak)

    # Information specific to broadenrich
	if(method == 'broadenrich' || method == 'broadenrich_splineless') {
		r_go_genes_avg_coverage = mean(gpw$ratio[b_genes])
	}

	# Small correction for case where every gene in this geneset has a peak.
	if(method == 'chipenrich_slow') {
		if (all(as.logical(sg_go))) {
			cont_length = quantile(gpw$length, 0.0025)

			cont_gene = data.frame(
				gene_id = "continuity_correction",
				length = cont_length,
				log10_length = log10(cont_length),
				num_peaks = 0,
				peak = 0,
				stringsAsFactors = FALSE)

			if ("mappa" %in% names(gpw)) {
				cont_gene$mappa = 1
			}
			gpw = rbind(gpw,cont_gene)
			b_genes = c(b_genes,1)

			message(sprintf("Applying correction for geneset %s with %i genes...", go_id, length(go_genes)))
		}
	}

    # Logistic regression works no matter the method because final_model is chosen above
    # and the data required from gpw will automatically be correct based on the method used.
    fit = gam(final_model, data = cbind(gpw, goterm = as.numeric(b_genes)), family = "binomial")

	# Results from the logistic regression
    r_effect = coef(fit)[2]
    r_pval = summary(fit)$p.table[2, 4]

	# The only difference between chipenrich and broadenrich here is
	# the Geneset Avg Gene Coverage column
	if(method == 'chipenrich_slow') {
		out = data.frame(
			"P.value" = r_pval,
			"Geneset ID" = r_go_id,
			"N Geneset Genes" = r_go_genes_num,
			"Geneset Peak Genes" = r_go_genes_peak,
			"N Geneset Peak Genes" = r_go_genes_peak_num,
			"Effect" = r_effect,
			"Odds.Ratio" = exp(r_effect),
			"Geneset Avg Gene Length" = r_go_genes_avg_length,
			stringsAsFactors = FALSE)
	} else if (method == 'broadenrich' || method == 'broadenrich_splineless') {
		out = data.frame(
			"P.value" = r_pval,
			"Geneset ID" = r_go_id,
			"N Geneset Genes" = r_go_genes_num,
			"Geneset Peak Genes" = r_go_genes_peak,
			"N Geneset Peak Genes" = r_go_genes_peak_num,
			"Effect" = r_effect,
			"Odds.Ratio" = exp(r_effect),
			"Geneset Avg Gene Length" = r_go_genes_avg_length,
			"Geneset Avg Gene Coverage" = r_go_genes_avg_coverage,
			stringsAsFactors = FALSE)
	}

	return(out)
}

test_gam = function(geneset,gpw,n_cores) {
	# Restrict our genes/weights/peaks to only those genes in the genesets.
	# Here, geneset is not all combined, but GOBP, GOCC, etc.
	# i.e. A specific one.
	gpw = subset(gpw, gpw$gene_id %in% geneset@all.genes)

	# Construct model formula.
	model = "peak ~ goterm + s(log10_length,bs='cr')"

	# Run tests. NOTE: If os == 'Windows', n_cores is reset to 1 for this to work
	results_list = parallel::mclapply(as.list(ls(geneset@set.gene)), function(go_id) {
		single_gam(go_id, geneset, gpw, 'chipenrich_slow', model)
	}, mc.cores = n_cores)

	# Collapse results into one table
	results = Reduce(rbind,results_list)

	# Correct for multiple testing
	results$FDR = stats::p.adjust(results$P.value, method = "BH")

	# Create enriched/depleted status column
	results$Status = ifelse(results$Effect > 0, 'enriched', 'depleted')

	results = results[order(results$P.value), ]

	return(results)
}

test_gam_ratio = function(geneset,gpw,n_cores) {
	# Restrict our genes/weights/peaks to only those genes in the genesets.
	gpw = subset(gpw, gpw$gene_id %in% geneset@all.genes)

	# Construct model formula.
	model = "goterm ~ ratio + s(log10_length,bs='cr')"

	# Run tests. NOTE: If os == 'Windows', n_cores is reset to 1 for this to work
	results_list = parallel::mclapply(as.list(ls(geneset@set.gene)), function(go_id) {
		single_gam(go_id, geneset, gpw, 'broadenrich', model)
	}, mc.cores = n_cores)

	# Collapse results into one table
	results = Reduce(rbind,results_list)

	# Correct for multiple testing
	results$FDR = stats::p.adjust(results$P.value, method = "BH")

	# Create enriched/depleted status column
	results$Status = ifelse(results$Effect > 0, 'enriched', 'depleted')

	results = results[order(results$P.value), ]

	return(results)
}

test_gam_ratio_splineless = function(geneset,gpw,n_cores) {
	message('Using test_gam_ratio_splineless..')

	# Restrict our genes/weights/peaks to only those genes in the genesets.
	gpw = subset(gpw, gpw$gene_id %in% geneset@all.genes)

	# Construct model formula.
	model = "goterm ~ ratio"

	# Run tests. NOTE: If os == 'Windows', n_cores is reset to 1 for this to work
	results_list = mclapply(as.list(ls(geneset@set.gene)), function(go_id) {
		single_gam(go_id, geneset, gpw, 'broadenrich_splineless', model)
	}, mc.cores = n_cores)

	# Collapse results into one table
	results = Reduce(rbind,results_list)

	# Correct for multiple testing
	results$FDR = stats::p.adjust(results$P.value, method = "BH")

	# Create enriched/depleted status column
	results$Status = ifelse(results$Effect > 0, 'enriched', 'depleted')

	results = results[order(results$P.value), ]

	return(results)
}
