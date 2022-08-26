#' faster version of 'test_chipenrich'
#' @keywords internal
#' @import methods
#' @import ChIPseeker
#' @import clusterProfiler
#' @importFrom IRanges distanceToNearest
#' @importFrom IRanges promoters
#' @author Brandon Signal
test_chipenrich_fast = function(geneset, gpw, n_cores=1, gpw.bg=NULL,bg_peaks=NULL, pre_filter=TRUE, genome, peak) {
  # Restrict our genes/weights/peaks to only those genes in the genesets.
  # Here, geneset is not all combined, but GOBP, GOCC, etc.
  # geneset in BG
  if(!is.null(bg_peaks)){
    geneset@all.genes = geneset@all.genes[geneset@all.genes %in% bg_peaks$gene_id]
  }
  gpw = subset(gpw, gpw$gene_id %in% geneset@all.genes)

  #gpw.bg = subset(gpw.bg, gpw.bg$gene_id %in% geneset@all.genes)

  # Making the first spline
  fitspl = mgcv::gam(peak~s(log10_length,bs='cr'),data=gpw,family="binomial")
  gpw$spline = as.numeric(predict(fitspl, gpw, type="terms"))

  # Construct model formula.
  model = "peak ~ goterm + spline"

  # Collapse results into one table
  #results.single = Reduce(rbind,results_list)

  #### FILTER WITH DIFFERENT ENRICHMENT METHOD FIRST
  if(pre_filter){
    message("Prefiltering genesets with enrichGO")

    ont = geneset@type
    if(ont == "Gene Ontology Biological Process") ont = "BP"
    if(ont == "Gene Ontology Molecular Function") ont = "MF"
    if(ont == "Gene Ontology Cellular Component") ont = "CC"

    if(!ont %in% c("CC", "BP", "MF")){
      warning("Ontology not matching for enrichGO. Ont=", ont)
    }

    if(genome == "mm10"){
      txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
      orgdb = org.Mm.eg.db::org.Mm.eg.db
    }else if(genome == "hg38"){
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      orgdb = org.Hs.eg.db::org.Hs.eg.db
    }else{
      stop("genome: ", genome, " not currently supported for use with chipenrich_fast")
    }

    # find nearest promotor (to mimic the "nearest_tss" method)
    dtn = IRanges::distanceToNearest(peak, IRanges::promoters(txdb, upstream = 1000, downstream = 0), ignore.strand=T)
    # get promotr locs
    proms = IRanges::promoters(txdb, upstream = 1000, downstream = 0)[dtn@to]
    # convert prom to gene
    gene <- ChIPseeker::seq2gene(proms, tssRegion = c(-1000, 1000), flankDistance = 0, TxDb=txdb)

    ego <- enrichGO(gene          = gene,
                    OrgDb         = orgdb,
                    ont           = ont,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01, #highish p-value (NOT ADJ!) cutoff
                    qvalueCutoff  = 0.5,
                    readable      = TRUE)
    ego=ego@result %>% filter(p.adjust < 0.05) #highish p-value adj cutoff

    test_go_ids = ego$ID
    if((length(test_go_ids)*4)>=length(ls(geneset@set.gene))){
      ego = arrange(ego, p.adjust)
      test_go_ids = ego$ID[1:(floor(length(ls(geneset@set.gene))/4))]
    }

    message("Filtered from ", length(ls(geneset@set.gene)), " to ", length(test_go_ids), " GO terms")

  }else{
    test_go_ids = NULL
  }

  results=NULL
  try({
    message("Using vectorised (chipenrichfast) enrichment function")
    results = multi_chipenrich(geneset, gpw, fitspl, method="chipenrich", model, skip_zero=T, test_go_ids = test_go_ids)
  })
  if(is.null(results)){
    message("Vectorisation failed... trying unvectorised...")
    # Run tests. NOTE: If os == 'Windows', n_cores is reset to 1 for this to work
    results_list = parallel::mclapply(as.list(test_go_ids), function(go_id) {
      single_chipenrich(go_id, geneset, gpw, 'chipenrich_slow', model, verbose = FALSE)
    }, mc.cores = n_cores)
    # Collapse results into one table
    results = Reduce(rbind,results_list)

    # Correct for multiple testing
    results$FDR = stats::p.adjust(results$P.value, method = "BH")


  }
  # Create enriched/depleted status column
  results$Status = ifelse(results$Effect > 0, 'enriched', 'depleted')

  results = results[order(results$P.value),]

  return(results)
}

#' faster version of 'single_chipenrich' that uses vectorised code'
#' @keywords internal
#' @import methods
#' @import ChIPseeker
#' @import clusterProfiler
#' @author Brandon Signal
multi_chipenrich = function(geneset, gpw, fitspl, method, model, skip_zero=T, test_go_ids=NULL, sample_size = 10, n_cores=1) {


  genes = ls(geneset@set.gene)
  if(!is.null(test_go_ids)){
    genes_filtered = genes[genes %in% test_go_ids]

    total_go_length = length(ls(geneset@set.gene))
    filtered_length = length(genes_filtered)
    missing_length = total_go_length - filtered_length

    #sample_size = 10
    repn = floor(missing_length / sample_size)

    genes_sampled = genes[!(genes %in% test_go_ids)]
    genes_sampled = sample(genes_sampled, sample_size)

    genes = c(genes_filtered, genes_sampled)
  }

  final_model = as.formula(model)


  # Genes in the geneset
  go_genes = lapply(genes, function(x) geneset@set.gene[[x]])
  # Background genes and the background presence of a peak
  b_genes = lapply(go_genes, function(x) gpw$gene_id %in% x)

  # Information about the geneset
  #r_go_id = go_id
  r_go_id = genes
  #r_go_genes_num = length(go_genes)
  r_go_genes_num = lapply(go_genes, length)


  # mean per GO geneset
  # SLOW AF
  #microbenchmark({
  #r_go_genes_avg_length = lapply(go_genes, function(x) mean(gpw$length[gpw$gene_id %in% x]))
  #}, unit = "ms", times=10L)

  go_df = data.frame(go_gene_id=unlist(go_genes), go_set=as.character(unlist(mapply(function(x, y) rep(x,length(y)), x=genes, y=go_genes))))
  go_df$gene_length = gpw$length[match(go_df$go_gene_id, gpw$gene_id)]
  r_go_genes_avg_length = aggregate(gene_length~go_set, go_df, mean)
  #r_go_genes_avg_length = mean(gpw$length[b_genes])

  # Information about peak genes
  go_df_peak = go_df[go_df$go_gene_id %in% unique(unlist(go_genes)) & go_df$go_gene_id %in% gpw$gene_id[gpw$peak==1],]
  r_go_genes_peak = aggregate(go_gene_id~go_set, go_df_peak, function(x) paste(x, collapse = ", "))
  r_go_genes_peak_num = aggregate(go_gene_id~go_set, go_df_peak, function(x) length(x))


  #go_genes_peak = gpw$gene_id[b_genes][sg_go==1]
  #r_go_genes_peak = paste(go_genes_peak,collapse=", ")
  #r_go_genes_peak_num = length(go_genes_peak)

  ## ONLY RUN IF ALL GENES HAVE A PEAK FOR GOTERM
  # Small correction for case where every gene in this geneset has a peak.

  go_df$peak = gpw$peak[match(go_df$go_gene_id, gpw$gene_id)]
  go_df_allpeaks = aggregate(peak ~ go_set, go_df, function(x) all(x==1))
  go_df_allpeaks = go_df_allpeaks[go_df_allpeaks$peak==TRUE,]

  ## FOR TESTING
  #go_df_allpeaks = aggregate(peak ~ go_set, go_df, function(x) length(which(x==1)))
  #go_df_allpeaks = go_df_allpeaks[which(go_df_allpeaks$peak ==113),]

  if (nrow(go_df_allpeaks) > 0) {
    cont_length = aggregate(gene_length~go_set, go_df[go_df$go_set %in% go_df_allpeaks$go_set,], function(x) quantile(x, 0.0025))
    #cont_length = quantile(gpw$length,0.0025)

    cont_gene = data.frame(
      gene_id = "continuity_correction",
      length = cont_length$gene_length,
      log10_length = log10(cont_length$gene_length),
      num_peaks = 0,
      peak = 0,
      go_set = cont_length$go_set,
      stringsAsFactors = FALSE)
    cont_gene$spline = as.numeric(predict(fitspl, cont_gene, type="terms"))

    if ("mappa" %in% names(gpw)) {
      cont_gene$mappa = 1
    }
    #gpw = rbind(gpw,cont_gene)
    #b_genes = c(b_genes,1)

    #message(sprintf("Applying correction for geneset %s with %i genes...",go_id,length(go_genes)))

    fix_peaks = which(genes %in% cont_length$go_set)
    b_genes_fix = b_genes[fix_peaks]
    b_genes = b_genes[-fix_peaks]

    out.fix = NULL
    for(i in 1:length(b_genes_fix)){

      try({
        #genes[fix_peaks][i]
        #cont_gene[cont_gene$go_set == genes[fix_peaks][i],-6]

        fit = mgcv::gam(final_model,data=cbind(rbind(gpw, cont_gene[cont_gene$go_set == genes[fix_peaks][i],-6]), goterm=c(as.numeric(b_genes_fix[[i]]), 1)),family="binomial")
        r_effect = coef(fit)[2];
        r_pval = summary(fit)$p.table[2, 4]

        out.fix = rbind(out.fix, data.frame(
          "P.value"=r_pval,
          "Effect"=r_effect,
          "Odds.Ratio"=exp(r_effect),
          stringsAsFactors = FALSE))
      })
    }
  }



  out.list = parallel::mclapply(b_genes, function(x){
    tryCatch({
      fit = mgcv::gam(final_model,data=cbind(gpw, goterm=as.numeric(x)),family="binomial")
      r_effect = coef(fit)[2];
      r_pval = summary(fit)$p.table[2, 4]

      out = data.frame(
        "P.value"=r_pval,
        "Effect"=r_effect,
        "Odds.Ratio"=exp(r_effect),
        stringsAsFactors = FALSE)
      return(out)
    },
    error = {function(e) {warning(
      #sprintf("Error in geneset: %s. NAs given", go_id))
      sprintf("Error in geneset: %s. NAs given", go_id))
    }}
    )

  }, mc.cores = n_cores)
  out = do.call("rbind", out.list)


  out2 = data.frame(
    "P.value"=out$P.value,
    "Geneset ID"=r_go_id,
    "N Geneset Genes"=unlist(r_go_genes_num),
    "Geneset Peak Genes"=r_go_genes_peak$go_gene_id[match(r_go_id, r_go_genes_peak$go_set)],
    "N Geneset Peak Genes"=r_go_genes_peak_num$go_gene_id[match(r_go_id, r_go_genes_peak_num$go_set)],
    "Effect"=out$Effect,
    "Odds.Ratio"=out$Odds.Ratio,
    "Geneset Avg Gene Length"=r_go_genes_avg_length$gene_length[match(r_go_id, r_go_genes_avg_length$go_set)],
    stringsAsFactors = FALSE)

  # fix for the genesets with 0 peaks
  out2$N.Geneset.Peak.Genes[which(is.na(out2$Geneset.Peak.Genes))] = 0
  out2$Geneset.Peak.Genes[which(is.na(out2$Geneset.Peak.Genes))] = ""

  if(!is.null(test_go_ids)){
    out2$source = ifelse(out2$Geneset.ID %in% genes_sampled, "SAMPLED", "FILTERED")
    out3 = rbind(out2, out2[rep(which(out2$source == "SAMPLED" & out2$P.value < 0.5), repn-1),])
    out3$FDR = p.adjust(out3$P.value, method="BH")
    out2 = out3[out3$source == "FILTERED",]
    out2$source = NULL
  }else{
    out2$FDR = p.adjust(out2$P.value, method="BH")
  }

  return(out2)
}

#' edited version of 'single_chipenrich' that can run quietly'
#' @keywords internal
#' @import methods
#' @author Brandon Signal

single_chipenrich = function(go_id, geneset, gpw, fitspl, method, model, verbose=T) {
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

  # Small correction for case where every gene in this geneset has a peak.
  if (all(as.logical(sg_go))) {
    cont_length = quantile(gpw$length,0.0025)

    cont_gene = data.frame(
      gene_id = "continuity_correction",
      length = cont_length,
      log10_length = log10(cont_length),
      num_peaks = 0,
      peak = 0,
      stringsAsFactors = FALSE)
    cont_gene$spline = as.numeric(predict(fitspl, cont_gene, type="terms"))


    if ("mappa" %in% names(gpw)) {
      cont_gene$mappa = 1
    }
    gpw = rbind(gpw,cont_gene)
    b_genes = c(b_genes,1)

    if(verbose) message(sprintf("Applying correction for geneset %s with %i genes...",go_id,length(go_genes)))
  }

  r_effect = NA
  r_pval = NA

  tryCatch(
    {fit = mgcv::gam(final_model,data=cbind(gpw,goterm=as.numeric(b_genes)),family="binomial")
    # Results from the logistic regression
    r_effect = coef(fit)[2];
    r_pval = summary(fit)$p.table[2, 4]
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
