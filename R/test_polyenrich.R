test_polyenrich = function(geneset, gpw, n_cores) {
	# Restrict our genes/weights/peaks to only those genes in the genesets.
	# Here, geneset is not all combined, but GOBP, GOCC, etc.
	# i.e. A specific one.
	gpw = subset(gpw, gpw$gene_id %in% geneset@all.genes)

	fitspl = mgcv::gam(num_peaks~s(log10_length,bs='cr'),data=gpw,family="nb")
	gpw$spline = as.numeric(predict(fitspl, gpw, type="terms"))

	# Construct model formula.
	model = "num_peaks ~ goterm + spline"

    nullfit = mgcv::gam(num_peaks~spline,data=gpw, family= "nb")

	# Run tests. NOTE: If os == 'Windows', n_cores is reset to 1 for this to work
	results_list = parallel::mclapply(as.list(ls(geneset@set.gene)), function(go_id) {
		single_polyenrich(go_id, geneset, gpw, fitspl, 'polyenrich', model, nullLR = mgcv::logLik.gam(nullfit))
	}, mc.cores = n_cores)

	# Collapse results into one table
	results = Reduce(rbind,results_list)

	# Correct for multiple testing
	results$FDR = p.adjust(results$P.value, method="BH")

	# Create enriched/depleted status column
	results$Status = ifelse(results$Effect > 0, 'enriched', 'depleted')

	results = results[order(results$P.value),]

	return(results)
}

single_polyenrich = function(go_id, geneset, gpw, fitspl, method, model, nullLR) {
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
	go_genes_peak = gpw$gene_id[b_genes][sg_go==1]
	r_go_genes_peak = paste(go_genes_peak,collapse=", ")
	r_go_genes_peak_num = length(go_genes_peak)

    r_effect = NA
    r_pval = NA

    tryCatch(
    {fit = mgcv::gam(final_model,data=cbind(gpw,goterm=as.numeric(b_genes)),family="nb")
        # Results from the logistic regression
        r_effect = coef(fit)[2];
        r_pval = stats::pchisq(2*(mgcv::logLik.gam(fit)-nullLR),1,lower.tail = F)
    },
    error = {function(e) {warning(
        sprintf("Error in geneset: %s. NAs given", go_id))
    }}
    )

	out = data.frame(
		"P.value"=r_pval,
		"Geneset ID"=r_go_id,
		"N Geneset Genes"=r_go_genes_num,
		"Geneset Peak Genes"=r_go_genes_peak,
		"N Geneset Peak Genes"=r_go_genes_peak_num,
		"Effect"=r_effect,
		"Odds.Ratio"=exp(r_effect),
		"Geneset Avg Gene Length"=r_go_genes_avg_length,
		stringsAsFactors = FALSE)

	return(out)
}
