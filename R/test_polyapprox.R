test_polyapprox = function(geneset, gpw, n_cores) {
	# Restrict our genes/weights/peaks to only those genes in the genesets.
	# Here, geneset is not all combined, but GOBP, GOCC, etc.
	# i.e. A specific one.
	gpw = subset(gpw, gpw$gene_id %in% geneset@all.genes)
	
    fitspl = mgcv::gam(num_peaks~s(log10_length,bs='cr'),data=gpw,family="nb")
    gpw$spline = as.numeric(predict(fitspl, gpw, type="terms"))
    
    #Need to use glm.nb to be able to use glm.scoretest later
    nullfit = MASS::glm.nb(num_peaks~spline,data=gpw)

	# Model formula not needed
	#model = "num_peaks ~ goterm + spline"
	
	# Run tests. NOTE: If os == 'Windows', n_cores is reset to 1 for this to work
	results_list = parallel::mclapply(as.list(ls(geneset@set.gene)), function(go_id) {
		single_polyapprox(go_id, geneset, gpw, 'polyenrich', nullmodel = nullfit)
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

glm.scoretest <- function(fit, x2, dispersion=NULL) #Directly lifted from statmod package,
#    because asking for more dependencies is asking for more chances for errors.
#    Score test for new covariate in glm
#    Gordon Smyth
#    27 March 2009. Last modified 20 Mar 2010.
{
    w <- fit$weights
    r <- fit$residuals
    if(any(w <= 0)) {
        r <- r[w>0]
        x2 <- x2[w>0]
        w <- w[w>0]
    }
    if (is.null(dispersion)) {
        fixed.dispersion <- (fit$family$family %in% c("poisson","binomial"))
        if(fixed.dispersion)
        dispersion <- 1
        else if(fit$df.residual > 0) {
            dispersion <- sum(w*r^2)/fit$df.residual
        } else {
            stop("No residual df available to estimate dispersion")
        }
    }
    ws <- sqrt(w)
    x2.1w <- qr.resid(fit$qr,ws*x2)
    zw <- ws*r
    colSums(as.matrix(x2.1w*zw))/sqrt(colSums(as.matrix(x2.1w * x2.1w)))/sqrt(dispersion)
}

single_polyapprox = function(go_id, geneset, gpw, method, nullmodel) {
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
    {r_effect = glm.scoretest(nullmodel, as.numeric(b_genes));
        r_pval = 2*pnorm(abs(r_effect),lower.tail = F)
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
