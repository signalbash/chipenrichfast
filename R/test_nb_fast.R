test_gam_nb_fast = function(geneset,gpw,n_cores) {
	# Restrict our genes/weights/peaks to only those genes in the genesets.
	# Here, geneset is not all combined, but GOBP, GOCC, etc.
	# i.e. A specific one.
	gpw = subset(gpw,gene_id %in% geneset@all.genes);

	if (sum(gpw$peak) == 0) {
		stop("Error: no peaks in your data!");
	}
	fitspl = gam(num_peaks~s(log10_length,bs='cr'),data=gpw,family="nb")
	gpw$spline = as.numeric(predict(fitspl, gpw, type="terms"))


	# Construct model formula.
	model = "num_peaks ~ goterm + spline";

	# Run tests on genesets and beware of Windows!
	if(os != 'Windows' && n_cores > 1) {
		results_list = mclapply(as.list(ls(geneset@set.gene)), function(go_id) {
			single_gam_nb_fast(go_id, geneset, gpw, fitspl, 'nbspeed', model)
		}, mc.cores = n_cores)
	} else {
		results_list = lapply(as.list(ls(geneset@set.gene)), function(go_id) {
			single_gam_nb_fast(go_id, geneset, gpw, fitspl, 'nbspeed', model)
		})
	}

	# Collapse results into one table
	results = Reduce(rbind,results_list)

	# Correct for multiple testing
	results$FDR = p.adjust(results$P.value, method="BH");

	# Create enriched/depleted status column
	results$Status = ifelse(results$Effect > 0, 'enriched', 'depleted')

	results = results[order(results$P.value),];

	return(results);
}

single_gam_nb_fast = function(go_id, geneset, gpw, fitspl, method, model) {
	final_model = as.formula(model);

	# Genes in the geneset
	go_genes = geneset@set.gene[[go_id]];

	# Filter genes in the geneset to only those in the gpw table.
	# The gpw table will be truncated depending on which geneset type we're in.
	go_genes = go_genes[go_genes %in% gpw$gene_id];

	# Background genes and the background presence of a peak
	b_genes = gpw$gene_id %in% go_genes;
	sg_go = gpw$peak[b_genes];

	# Information about the geneset
	r_go_id = go_id;
	r_go_genes_num = length(go_genes);
	r_go_genes_avg_length = mean(gpw$length[b_genes]);

	# Information about peak genes
	go_genes_peak = gpw$gene_id[b_genes][sg_go==1];
	r_go_genes_peak = paste(go_genes_peak,collapse=", ");
	r_go_genes_peak_num = length(go_genes_peak);

	# Logistic regression works no matter the method because final_model is chosen above
	# and the data required from gpw will automatically be correct based on the method used.
	fit = gam(final_model,data=cbind(gpw,goterm=as.numeric(b_genes)),family="nb");

	# Results from the logistic regression
	r_effect = coef(fit)[2];
	r_pval = summary(fit)$p.table[2,4];

	# The only difference between chipenrich and broadenrich here is
	# the Geneset Avg Gene Coverage column

	out = data.frame(
		"P.value"=r_pval,
		"Geneset ID"=r_go_id,
		"N Geneset Genes"=r_go_genes_num,
		"Geneset Peak Genes"=r_go_genes_peak,
		"N Geneset Peak Genes"=r_go_genes_peak_num,
		"Effect"=r_effect,
		"Odds.Ratio"=exp(r_effect),
		"Geneset Avg Gene Length"=r_go_genes_avg_length,
		stringsAsFactors = FALSE);

	return(out);
}
