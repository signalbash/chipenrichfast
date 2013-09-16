#!/usr/bin/env Rscript
library(chipenrich.data);
library(mgcv);
library(IRanges);
library(GenomicRanges);
library(lattice);
library(latticeExtra);
library(grid);
library(stringr);
library(reshape);

SUPPORTED_METHODS = list(
  'chipenrich' = "test_gam",
  'fet' = "test_fisher_exact"
);

HIDDEN_METHODS = list(
  'binomial' = "test_binomial"
);

METHOD_NAMES = list(
	'fet' = "Fisher's Exact Test",
	'binomial' = "Binomial Test",
	'chipenrich' = "ChIP-Enrich"
);

read_bed = function(file_path) {
  if (!file.exists(file_path)) {
    stop("Can't find BED file: ",file_path);
  }
  
  chunk_size = 500;
  chunk = scan(file_path,what="character",nmax=chunk_size,strip.white=T,sep="\n",quiet=T);
  skip_n = suppressWarnings(min(grep("^chr(\\d+|\\w+)\\s+\\d+",chunk)) - 1);
  
  if (is.infinite(skip_n)) {
    stop("Error: no valid chromosomes detected within first 500 lines of BED file.");
  }
  
  message(sprintf("Skipping %i lines of BED header..",skip_n));
  
  peaks = read.table(file_path,header=F,skip=skip_n);
  peaks = peaks[,1:3];
  names(peaks) = c("chrom","start","end");

  sub_check = peaks[1:min(nrow(peaks),100),];
  if (!all(grepl("chr",sub_check$chrom))) {
    stop("First column of BED file should have chr* entries.");
  }
  
  if (any(sub_check$start < 0) | any(sub_check$end < 0)) {
    stop("Start/end positions of peaks should be >= 0.");
  }
  
  chroms = list();
  for (chr in unique(peaks$chrom)) {
    peaks_chrom = subset(peaks,chrom == chr);
    chroms[[chr]] = IRanges(start=peaks_chrom$start,end=peaks_chrom$end);
  }
  
  return(chroms);
}

read_peaks = function(file_path) { 
  d = read.table(file_path,sep="\t",header=T);
  
  # Check columns. 
  for (col in c("chrom","start","end")) {
    if (!col %in% names(d)) {
      stop(sprintf("error reading peaks: no column named '%s' in file",col));
    }
  }
  
  chroms = list();
  for (chr in unique(d$chrom)) {
    peaks_chrom = subset(d,chrom == chr);
    chroms[[chr]] = IRanges(start=peaks_chrom$start,end=peaks_chrom$end);
  }
  
  return(chroms);
}

load_peaks = function(dframe) {
  # Check columns. 
  for (col in c("chrom","start","end")) {
    if (!col %in% names(dframe)) {
      stop(sprintf("error in peaks data frame: no column named '%s'",col));
    }
  }
  
  chroms = list();
  for (chr in unique(dframe$chrom)) {
    peaks_chrom = subset(dframe,chrom == chr);
    chroms[[chr]] = IRanges(start=peaks_chrom$start,end=peaks_chrom$end);
  }
  
  return(chroms);
}

assign_peaks = function(peaks,locusdef,tss,midpoint=T) {
  results = list();
  
  if (!midpoint) {
    warning("peak midpoints must be used currently - may be changed in future versions");
  }
  
  for (chrom in names(peaks)) {
    if (!chrom %in% names(locusdef@chrom2iranges)) {
      next;
    }
    
    peak_ranges = peaks[[chrom]];
    
    # Use midpoints of peak ranges. 
    peak_mids = mid(peak_ranges);
    peak_mid_ranges = IRanges(peak_mids,peak_mids);
    
		# Find overlapping peaks. 
    gene_ranges = locusdef@chrom2iranges[[chrom]];
    overlaps = findOverlaps(peak_mid_ranges,gene_ranges);
    overlap_matrix = IRanges::as.matrix(overlaps);

    if (length(overlaps) == 0) {
			# There were no peak overlaps with gene loci on this chromosome. 
      next;
    }
    
    matched_orig = peak_ranges[overlap_matrix[,1]];
    matched_peaks = peak_mid_ranges[overlap_matrix[,1]]; # Peaks that had an overlap with a gene locus
    matched_genes = gene_ranges[overlap_matrix[,2]]; # The gene locus assigned to each peak
    
		# Find the nearest TSS and the distance to it from each matched peak (above). 
		# tss_ranges = tss$ranges[[chrom]];
    # nearest_tss = tss_ranges[nearest(matched_peaks,tss_ranges)];
    # dist_to_tss = distance(matched_peaks,nearest_tss);
		
    # Find nearest TSS, distance to it (signed by up (-) or down (+) stream. 
    matched_peak_grange = GRanges(chrom,matched_peaks);
    tss_hits = tss$granges[nearest(matched_peak_grange,tss$granges)];
    
    neg1 = start(matched_peak_grange) < start(tss_hits) & (strand(tss_hits) == "+");
    neg2 = start(matched_peak_grange) > start(tss_hits) & (strand(tss_hits) == "-");
    neg_dist = as.logical(neg1 | neg2);
    
    suppressWarnings({
      # Hide the Bioconductor 2.12 warning about distance functionality changing
      dist_obj = distanceToNearest(matched_peak_grange,tss$granges);
    });

    if (class(dist_obj) == "Hits") {
      dist_to_tss = dist_obj@elementMetadata$distance;
    } else {
      dist_to_tss = dist_obj$distance;
    }
    
    dist_to_tss[neg_dist] = dist_to_tss[neg_dist] * -1;
    
    d = data.frame(
      chrom=chrom,
      peak_start=start(matched_orig),
      peak_end=end(matched_orig),
      peak_midpoint=start(matched_peaks),
      gene_locus_start=start(matched_genes),
      gene_locus_end=end(matched_genes),
      geneid=names(matched_genes),
      nearest_tss=start(tss_hits),
      nearest_tss_gene=tss_hits$geneid,
      dist_to_tss=dist_to_tss,
      nearest_tss_gene_strand=as.character(strand(tss_hits))
    );
    
		# if (midpoint) {
			# d = data.frame(
				# chrom=chrom,
				# peak_midpoint=start(matched_peaks),
				# locus_start=start(matched_genes),
				# locus_end=end(matched_genes),
				# geneid=names(matched_genes),
        # nearest_tss=start(tss_hits),
        # dist_to_tss=dist_to_tss$distance,
        # strand=as.character(strand(tss_hits))
			# );
		# } else {
			# d = data.frame(
				# chrom=chrom,
				# peak_start=start(matched_peaks),
				# peak_end=end(matched_peaks),
				# locus_start=start(matched_genes),
				# locus_end=end(matched_genes),
				# dist_to_tss=dist_to_tss$distance,
				# geneid=names(matched_genes)
			# );
		# }
    
    results[[chrom]] = d;
  }
  
  if (length(results) == 0) {
    return(NULL);
  } else{
    result = merge_all(results);
    return(result);
  }
}

num_peaks_per_gene = function(assigned_peaks,locusdef,mappa=NULL) { 
  # Add in gene lengths to locusdef data. 
  d = locusdef@dframe;
  d$length = d$end - d$start;
  
  # Sum up lengths for each gene. 
  d = stats::aggregate(length ~ geneid,d,sum);
  d$log10_length = log10(d$length);
  
  # Compute the total number of peaks assigned to each gene. 
  ppg = table(assigned_peaks$geneid);
  d_ppg = data.frame(geneid=names(ppg),num_peaks=as.numeric(ppg),stringsAsFactors=F);
  result = merge(d,d_ppg,by="geneid",all.x=T);
  result[is.na(result$num_peaks),]$num_peaks = 0;
  
  # Mappable length if requested. 
	if (!is.null(mappa)) {
    message("Mappability adjustment is enabled..");
    result = merge(result,mappa,by="geneid",sort=F);
    result$length = as.numeric((result$mappa * result$length) + 1);
    result$log10_length = log10(result$length);
  }
  
  # Order by number of peaks in a gene. 
  result = result[order(result$num_peaks,decreasing=T),];
  
  # Add in peak vector (0,1). 
  result$peak = as.numeric(result$num_peaks >= 1);
  
  return(result);
}

add_emp_peak = function(gpw,bin_size=25) {
  d = gpw;
  d = d[order(d$log10_length),];
  d$group = ceiling((1:dim(d)[1])/bin_size);
  
  bygroup = stats::aggregate(peak ~ group,d,mean);
  d$emp_peak = bygroup$peak[d$group];
  d;
}

add_emp_peak_mappa = function(gpw,bin_size=25) {
  d = gpw;
  d = d[order(d$mappa),];
  d$group = ceiling((1:dim(d)[1])/bin_size);
  
  bygroup = stats::aggregate(peak ~ group,d,mean);
  d$emp_peak = bygroup$peak[d$group];
  d;
}

avg_binned_peak = function(gpw,bin_size=25) {
  d = gpw;
  d = d[order(d$log10_length),];
  d$group = ceiling((1:dim(d)[1])/bin_size);
  
  bygroup = stats::aggregate(cbind(peak,length) ~ group,d,mean);
  bygroup$log_avg_length = log10(bygroup$length);
  names(bygroup) = c("group","peak","avg_length","log_avg_length");
  bygroup;
}

avg_binned_peak_mappa = function(gpw,bin_size=25) {
  d = gpw;
  d = d[order(d$mappa),];
  d$group = ceiling((1:dim(d)[1])/bin_size);
  
  bygroup = stats::aggregate(cbind(peak,mappa) ~ group,d,mean);
  bygroup$log_avg_mappa = log10(bygroup$mappa);
  names(bygroup) = c("group","peak","avg_mappa","log_avg_mappa");
  bygroup;
}

make_gpw = function(locusdef,peak_genes) {
  d = locusdef@dframe;
  
  # Indicator vector for which genes have peaks. 
  d$peak = as.numeric(d$geneid %in% peak_genes);
  
  # Compute length and log10 length for each gene. 
  d$length = d$end - d$start;
  
  # For genes that exist across multiple rows, sum their lengths. 
  d = stats::aggregate(cbind(peak,length) ~ geneid,d,sum);
  d$log10_length = log10(d$length);
  
  # A gene could now have > 1 peak due to aggregating (above), reset it to be 
  # 1 or 0. 
  d$peak = as.integer(d$peak >= 1);

  # Sort by locus length. 
  d = d[order(d$log10_length),];

  return(subset(d,select=c("geneid","length","log10_length","peak")));
}

calc_weights_gam = function(locusdef,peak_genes,mappa=NULL,...) {
  d = locusdef@dframe;
  
  # Indicator vector for which genes have peaks. 
  d$peak = as.numeric(d$geneid %in% peak_genes);
  
  # Compute length and log10 length for each gene. 
  d$length = d$end - d$start;
  
  # For genes that exist across multiple rows, sum their lengths. 
  d = stats::aggregate(cbind(peak,length) ~ geneid,d,sum);
  d$log10_length = log10(d$length);
  
  # A gene could now have > 1 peak due to aggregating (above), reset it to be 
  # 1 or 0. 
  d$peak = as.integer(d$peak >= 1);
  
  # Sort by locus length. 
  d = d[order(d$log10_length),];
  
  # If mappability was requested, add that in. 
  if (!is.null(mappa)) {
    d = merge(d,mappa,by="geneid",sort=F);
    d$orig_length = d$length;
    d$length = as.numeric((d$mappa * d$length) + 1);
    d$log10_length = log10(d$length);
    d = d[order(d$log10_length),];
  }
  
  # Create model. 
  model = "peak ~ s(log10_length,bs='cr')";
  
  # Compute binomial spline fit.
  fit = gam(as.formula(model),data=d,family="binomial");
  
  # Compute weights for each gene, based on the predicted prob(peak) for each gene. 
  ppeak = fitted(fit);
  w0 = 1 / (ppeak/mean(d$peak,na.rm=T));
  w0 = w0 / mean(w0,na.rm=T);
  
  d$weight = w0;
  d$prob_peak = ppeak;
  d$resid.dev = resid(fit,type="deviance");
  
  cols = c("geneid","length","log10_length","mappa","orig_length","peak","weight","prob_peak","resid.dev");
  if (is.null(mappa)) {
    cols = setdiff(cols,c("mappa","orig_length"));
  }
  d = subset(d,select=cols);
  return(d);
}

calc_spline_mappa = function(locusdef,peak_genes,mappa,...) {
  d = locusdef@dframe;
  
  # Indicator vector for which genes have peaks. 
  d$peak = as.numeric(d$geneid %in% peak_genes);
  
  # Compute length and log10 length for each gene. 
  d$length = d$end - d$start;
  
  # For genes that exist across multiple rows, sum their lengths. 
  d = stats::aggregate(cbind(peak,length) ~ geneid,d,sum);
  d$log10_length = log10(d$length);
  
  # A gene could now have > 1 peak due to aggregating (above), reset it to be 
  # 1 or 0. 
  d$peak = as.integer(d$peak >= 1);
  
  # Merge in mappability.  
  d = merge(d,mappa,by="geneid",sort=F);
  
  # Sort by mappa. 
  d = d[order(d$mappa),];
  
  # Create model. 
  model = "peak ~ s(mappa,bs='cr')";
  
  # Compute binomial spline fit.
  fit = gam(as.formula(model),data=d,family="binomial");
  
  # Compute weights for each gene, based on the predicted prob(peak) for each gene. 
  d$prob_peak = fitted(fit);
  return(subset(d,select=c("geneid","length","log10_length","mappa","peak","prob_peak")));
}

test_fisher_exact = function(geneset,gpw,alternative="two.sided") {
  # Restrict to only those genes in the genesets. 
  gpw = subset(gpw,geneid %in% geneset@all.genes);

  genes = gpw$geneid;
  peaks = gpw$peak;

  results = lapply(ls(geneset@set.gene), function(go_term) {
    go_genes = geneset@set.gene[[go_term]];
    n_go_genes = length(go_genes);

    in_cat = as.numeric(genes %in% go_genes);
    n_go_peak_genes = sum((genes %in% go_genes) & (peaks == 1));
    xt = table(in_cat,peaks);    
    
    pval = 1;
    odds_ratio = 0;
    try({ 
      fet_result = fisher.test(xt,conf.int=F,alternative=alternative);
      pval = fet_result$p.value;
      odds_ratio = fet_result$estimate;
    },silent=T);
    
		go_peak_genes = gpw[(gpw$geneid %in% go_genes) & (gpw$peak == 1),]$geneid;
		go_peak_genes_str = paste(go_peak_genes,collapse=", ");
		
    enr = NA;
    
    if (odds_ratio > 1) {
      enr = "enriched";
    }
    
    if (odds_ratio < 1) {
      enr = "depleted";
    }
    
    if (odds_ratio == 0) { 
      enr = "depleted";
    }
    
    if (is.infinite(odds_ratio)) {
      enr = "enriched";
    }
    
    data.frame(
      "Geneset ID" = go_term,
      "N Geneset Genes" = n_go_genes,
      "N Geneset Peak Genes" = n_go_peak_genes,
			"Geneset Peak Genes" = go_peak_genes_str,
      "Odds.Ratio" = odds_ratio,
      "Status" = enr,
      "P-value" = pval,
      stringsAsFactors=F
    );
  });

  results = rbind.fill(results);  

  results$FDR = p.adjust(results$P.value,method="BH");

  results = results[order(results$P.value),];
  return(results);
}

test_binomial = function(geneset,ppg) {
  # Restrict to genes only within the genesets. 
  ppg = subset(ppg,geneid %in% geneset@all.genes);
  
  # Total # of peaks. 
  total_peaks = sum(ppg$num_peaks);

  # If there are no total peaks, it means there were no genes w/ peaks that overlap
  # with all genes across all gene sets. We can't analyze in this case. 
  if (total_peaks == 0) {
    stop("None of the peaks in this dataset were assigned to any of the genes across all pathways/genesets being tested.");
  }

  # Gene IDs must be character to be stored into an environment as names. 
  ppg$geneid = as.character(ppg$geneid);

  # Index gene --> peak count. 
  gene2peakc = new.env(parent=emptyenv());
  apply(ppg,1,function(x) assign(x["geneid"],x["num_peaks"],envir=gene2peakc));
  
  # Index gene --> log10 length. 
  gene2length = new.env(parent=emptyenv());
  apply(ppg,1,function(x) assign(x["geneid"],x["length"],envir=gene2length));

  # Length of the genome. 
  genome_length = sum(as.numeric(ppg$length));
  
  estimates = c();  
  pvals = c();
  goids = c();
  gonums = c();
  n_go_peak_genes = c();
	r_go_peak_genes = c();
  for (i in 1:length(geneset@set.gene)) {
    go_term = ls(geneset@set.gene)[i];
    go_genes = as.character(geneset@set.gene[[go_term]]);
		
		# Eliminate GO genes that aren't in the ppg. 
		#go_genes = Filter(function(x) x %in% ppg$geneid,go_genes);
    go_genes = go_genes[go_genes %in% ppg$geneid];

    go_gene_lengths = as.numeric(sapply(go_genes,function(x) get(x,gene2length)));
    go_term_length = sum(go_gene_lengths);

    go_gene_counts = as.numeric(sapply(go_genes,function(x) get(x,gene2peakc)));
    go_term_count = sum(go_gene_counts);    

    binom_p = go_term_length / genome_length;
    btest = binom.test(go_term_count,total_peaks,binom_p,alternative="greater");  
  
    goids[i] = go_term;
    pvals[i] = btest$p.value;
    estimates[i] = btest$estimate;
    gonums[i] = length(go_genes);
		
		go_peak_genes = ppg[(ppg$geneid %in% go_genes) & (ppg$num_peaks > 0),]$geneid;
		r_go_peak_genes[i] = paste(go_peak_genes,collapse=", ");
		n_go_peak_genes[i] = length(go_peak_genes);
  }
  
  fdr = p.adjust(pvals,method="BH");
  
  results = data.frame(
    "Geneset ID"=goids,
    "N Geneset Genes"=gonums,
    "N Geneset Peak Genes"=n_go_peak_genes,
    "P.Success"=estimates,
    "P-value"=pvals,
    "FDR"=fdr,
		"Geneset Peak Genes"=r_go_peak_genes,
    stringsAsFactors=F
  );
  
  results = results[order(results$P.value),];
  return(results);
}

test_gam = function(geneset,gpw) { 
  # Restrict our genes/weights/peaks to only those genes in the genesets. 
  gpw = subset(gpw,geneid %in% geneset@all.genes);

  if (sum(gpw$peak) == 0) {
    stop("Error: no peaks in your data!");
  }

  # Construct model formula. 
  model = "peak ~ goterm + s(log10_length,bs='cr')";
  final_model = as.formula(model);
    
  r_pvals = c();
  r_go_ids = c();
  r_go_genes = c();
  r_go_genes_num = c();
  r_go_genes_peak = c();
  r_go_genes_peak_num = c();
  r_effects = c();
  for (i in 1:length(geneset@set.gene)) {
    go_id = ls(geneset@set.gene)[i];
    go_genes = geneset@set.gene[[go_id]];
		
		# Eliminate GO genes that aren't in the ppg. 
		#go_genes = Filter(function(x) x %in% gpw$geneid,go_genes);
    go_genes = go_genes[go_genes %in% gpw$geneid];
    
    b_genes = gpw$geneid %in% go_genes;
    sg_go = gpw$peak[b_genes];

    if ((sum(b_genes) == 0) || (sum(!b_genes) == 0)) {
      next;
    }   
 
    r_go_genes[i] = paste(go_genes,collapse=";");
    r_go_genes_num[i] = length(go_genes);
    
    go_genes_peak = gpw$geneid[b_genes][sg_go==1];
    r_go_genes_peak[i] = paste(go_genes_peak,collapse=", ");
    r_go_genes_peak_num[i] = length(go_genes_peak);
    
    if (all(as.logical(sg_go))) {
      # Small correction for case where every gene in this geneset has a peak. 
      cont_length = quantile(gpw$length,0.0025);
      cont_gene = data.frame(
        geneid = "continuity_correction",
        length = cont_length,
        log10_length = log10(cont_length),
        num_peaks = 0,
        peak = 0,
        stringsAsFactors = F
      );
      if ("mappa" %in% names(gpw)) {
        cont_gene$mappa = 1;
      }
      gpw = rbind(gpw,cont_gene);
      b_genes = c(b_genes,1);

      message(sprintf("Applying correction for geneset %s with %i genes...",go_id,length(go_genes)));
    }
    
    fit = gam(final_model,data=cbind(gpw,goterm=as.numeric(b_genes)),family="binomial");
    effect = coef(fit)[2];
    bicep = summary(fit)$p.table[2,4];

    r_effects[i] = effect;
    r_pvals[i] = bicep;
    r_go_ids[i] = go_id;
  }
  
  fdr = p.adjust(r_pvals,method="BH");
  
  is_depleted = r_effects < 0;
  is_enriched = r_effects > 0;
  
  enr = rep(NA,length(r_go_ids));
  enr[is_depleted] = "depleted";
  enr[is_enriched] = "enriched";
  
  results = data.frame(
    "Geneset ID"=r_go_ids,
    "N Geneset Genes"=r_go_genes_num,
    "N Geneset Peak Genes"=r_go_genes_peak_num,
    "Effect"=r_effects,
    "Odds.Ratio"=exp(r_effects),
    "Status"=enr,
    "P-value"=r_pvals,
    "FDR"=fdr,
		"Geneset Peak Genes"=r_go_genes_peak,
    stringsAsFactors=F
  );
  
  results = results[order(results$P.value),];
  return(results);
}

# Quick function to change column names of data frame 
# (by name, not by position.) 
change_names = function(dframe,name_list) {
  for (n in names(name_list)) {
    names(dframe)[names(dframe)==n] = name_list[[n]];
  }
  dframe;
}

# Remove a file extension from a file name. 
strip_ext = function(file) {
  s = unlist(strsplit(file,".",fixed=T));
  s = s[1:(length(s)-1)];
  s = paste(s,collapse=".");
  return(s);
}

# Get a file's extension. 
get_ext = function(file) { 
  s = unlist(strsplit(file,".",fixed=T));
  return(tail(s,1));
}

..plot_expected_peaks = function(ppg) {
  # Gene/peak counts. 
  num_genes = dim(ppg)[1];
  num_peaks = sum(ppg$num_peaks);

  # Compute prob of mappa for each gene. 
  genome_length = sum(as.numeric(ppg$length));
  ppg$prob_mappable = ppg$length / genome_length;
  
  # Expected peaks (lambda)
  ppg$exp_peaks = num_peaks * ppg$prob_mappable;

  # Quantiles. 
  ppg$qpois_05 = qpois(0.05,ppg$exp_peaks);
  ppg$qpois_95 = qpois(0.95,ppg$exp_peaks);
  ppg$qpois_bonf_05 = qpois(0.05 / num_genes,ppg$exp_peaks);
  ppg$qpois_bonf_95 = qpois(1 - (0.05 / num_genes),ppg$exp_peaks);
  
  # Max # of peaks, expected or observed. 
  max_peaks = max(ppg$num_peaks,ppg$expected_peaks);
  pad = round(0.05 * max_peaks);
  
  # Thin out num_peaks and expected_peaks vs. length. 
  d_npeak = ppg[,c("num_peaks","length")];
  d_npeak = unique(d_npeak);
  
  d_epeak = ppg[,c("exp_peaks","length")];
  d_epeak = unique(d_epeak);
  
  # Colors. 
  col_npeaks = "blue";
  col_epeaks = "black";
  col_qpois_05 = "green";
  col_qpois_95 = "darkgreen";
  col_qpois_bonf_05 = "red";
  col_qpois_bonf_95 = "darkred";

  custom_key = list(
    text = list(
      c(
        'Count of Peaks',
        'Expected Peaks',
        '5th Percentile',
        '95th Percentile',
        '5th Percentile (Bonferroni)',
        '95th Percentile (Bonferroni)'
      )
    ),
    lines = list(
      pch = c(20),
      type = c("p","l","l","l","l","l"),
      lwd = 2,
      col = c(
        col_npeaks,
        col_epeaks,
        col_qpois_05,
        col_qpois_95,
        col_qpois_bonf_05,
        col_qpois_bonf_95
      )
    ),
    cex = 1.2,
    columns = 2
  );
  
  p = xyplot(
    num_peaks ~ length,
    d_npeak,
    par.settings = simpleTheme(pch=20,col=col_npeaks),
    ylim = c(-pad,max_peaks+pad),
    key = custom_key,
    scales = list(cex=1.4),
    xlab = list(label="Gene Locus Length",cex=1.4),
    ylab = list(label="Count of Peaks",cex=1.4)
  );
  
  p = p + as.layer(xyplot(
    exp_peaks ~ length,
    d_epeak,
    par.settings = simpleTheme(pch=20,col=col_epeaks,lwd=2),
    type = "l"
  ));
  
  for (v in c("qpois_05","qpois_95","qpois_bonf_05","qpois_bonf_95")) {
    d = ppg[,c(v,"length")];
    d = unique(round(d,2));
    d = d[order(d$length),];
    
    pcol = get(sprintf("col_%s",v));
    form = as.formula(sprintf("%s ~ length",v));
    p = p + as.layer(xyplot(
      form,
      d,
      par.settings = simpleTheme(pch=20,col=pcol,lwd=2),
      type = "l"
    ));
  }
  
  # Plot. 
  # p = xyplot(
    # num_peaks + exp_peaks + qpois_05 + qpois_95 + qpois_bonf_05 + qpois_bonf_95 ~ length,
    # ppg,
    # par.settings = simpleTheme(pch=20),
    # auto.key = list(
      # cex = 1.3,
      # text = c(
        # 'Count of Peaks',
        # 'Expected Peaks',
        # '5th Percentile',
        # '95th Percentile',
        # '5th Percentile (Bonferroni)',
        # '95th Percentile (Bonferroni)'
      # ),
      # columns = 2
    # ),
    # scales = list(cex=1.4),
    # xlab = list(label="Gene Length",cex=1.4),
    # ylab = list(label="Count of Peaks",cex=1.4)
  # );
  
  return(p);
}

plot_expected_peaks = function(peaks,locusdef="nearest_tss",genome="hg19",use_mappability=F,read_length=36,mappa_file=NULL) {
   # Check genome. 
  if (!genome %in% supported_genomes()) {
    stop("genome not supported: ",genome);
  }
  
  # Organism - human or mouse. 
  organism = genome_to_organism(genome);
  
  # Check locus definition. Should only be 1. 
  if (!locusdef %in% supported_locusdefs()) {
    stop("bad locus definition requested: ",locusdef);
  }
  
  # Check read length. 
  if (use_mappability) {
    if (!as.numeric(read_length) %in% supported_read_lengths()) {
      stop("bad read length requested: ",read_length);
    }
  }
  
  # Get peaks from user's file. 
	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks);
	} else if (class(peaks) == "character") {
    if (get_ext(peaks) == "bed") {
      message("Reading BED file: ",peaks);
      peakobj = read_bed(peaks);
    } else {
      message("Reading peaks file: ",peaks);
      peakobj = read_peaks(peaks);
    }
	}
  
	# Number of peaks in data. 
  num_peaks = sum(sapply(peakobj,function(x) length(x)))
  
  # Load locus definitions. 
  ldef_code = sprintf("locusdef.%s.%s",genome,locusdef);
  data(list=ldef_code,package = "chipenrich.data");
  ldef = get(ldef_code);
  
  # Load mappability if requested. 
  if (!is.null(mappa_file)) {
    message("Reading user-defined gene locus mappability file: ",mappa_file);
    mappa = read_mappa(mappa_file);
  } else {
    if (use_mappability) {
      mappa_code = sprintf("mappa.%s.%s.%imer",genome,locusdef,read_length);
      data(list=mappa_code,package = "chipenrich.data");
      mappa = get(mappa_code);
      mappa = na.omit(mappa);
    } else {
      mappa = NULL;
    }
  }
  
  # Load TSS site info. 
  tss_code = sprintf("tss.%s",genome);
  data(list=tss_code,package = "chipenrich.data");
  tss = get(tss_code);
  
  # Assign peaks to genes. 
  message("Assigning peaks to genes..");
  assigned_peaks = assign_peaks(peakobj,ldef,tss,midpoint=T);
  
  # Data about the # of peaks per gene. 
  ppg = num_peaks_per_gene(assigned_peaks,ldef,mappa);
  
  # Plot. 
  p = ..plot_expected_peaks(ppg);
  
  return(p);
}

plot_spline_length = function(peaks,locusdef="nearest_tss",genome='hg19',use_mappability=F,read_length=36,legend=T,xlim=NULL) {
	# Check genome. 
  if (!genome %in% supported_genomes()) {
    stop("genome not supported: ",genome);
  }
  
  # Check locus definition. Should only be 1. 
  if (!locusdef %in% supported_locusdefs()) {
    stop("bad locus definition requested: ",locusdef);
  }
  
  # Check read length. 
  if (use_mappability) {
    if (!as.numeric(read_length) %in% supported_read_lengths()) {
      stop("bad read length requested: ",read_length);
    }
  }
	
	# Get peaks from user's file. 
	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks);
	} else if (class(peaks) == "character") {
    if (get_ext(peaks) == "bed") {
      message("Reading BED file: ",peaks);
      peakobj = read_bed(peaks);
    } else {
      message("Reading peaks file: ",peaks);
      peakobj = read_peaks(peaks);
    }
	}
  
	# Number of peaks in data. 
  num_peaks = sum(sapply(peakobj,function(x) length(x)))
  
  # Load locus definitions. 
  ldef_code = sprintf("locusdef.%s.%s",genome,locusdef);
  data(list=ldef_code,package = "chipenrich.data");
  ldef = get(ldef_code);
  
	# Load TSS site info. 
  tss_code = sprintf("tss.%s",genome);
  data(list=tss_code,package = "chipenrich.data");
  tss = get(tss_code);
  
  # Load mappability if requested. 
  if (use_mappability) {
    mappa_code = sprintf("mappa.%s.%s.%imer",genome,locusdef,read_length);
    data(list=mappa_code,package = "chipenrich.data");
    mappa = get(mappa_code);
  } else {
    mappa = NULL;
  }
	
  # Assign peaks to genes. 
  assigned_peaks = assign_peaks(peakobj,ldef,tss);
  peak_genes = unique(assigned_peaks$geneid);
  
	# Make plot. 
	plotobj = ..plot_spline_length(ldef,peak_genes,num_peaks,mappa,legend=legend,xlim=xlim);
	return(plotobj);
}

# Create diagnostic plot of P(peak) from your data against log locus length (of genes)
# along with the spline fit.  
..plot_spline_length = function(locusdef,peak_genes,num_peaks,mappa=NULL,legend=T,xlim=NULL) {  
  # Calculate smoothing spline fit.
  gpw = calc_weights_gam(locusdef,peak_genes,mappa=mappa); # gpw = genes, peaks, weights
  
  # Genome length. 
  genome_length = sum(as.numeric(gpw$length));
    
  # Add in empirical peak by binning. 
  gpw = add_emp_peak(gpw);
  
  # Average peak/lengths. 
  avg_bins = avg_binned_peak(gpw,bin_size=25);
  
  # Order by length. 
  gpw = gpw[order(gpw$log10_length),];
  
  # Calculate prob(all false positives). 
  gpw$false_prob = 1 - (1 - (gpw$length/genome_length))^(num_peaks);
  
  col_rand_gene = "grey35";
  col_rand_peak = "grey74";
  col_spline = "darkorange";
  
  panel_func = function(x,y,...) { 
    panel.xyplot(x,y,...); 
    panel.abline(h=mean(y),col=col_rand_gene,lwd=3);
  }

  custom_key = list( 
    text = list(
      c(
        "Expected Fit - Peaks Independent of Locus Length",
        "Expected Fit - Peaks Proportional to Locus Length",
        "Binomial Smoothing Spline Fit",
        "Proportion of Genes in Bin with at Least 1 Peak"
      )
    ),
    lines = list(
      pch = c(20,15,17,20),
      type = c("l","l","l","p"),
      col = c(col_rand_gene,col_rand_peak,col_spline,"black"),
      lwd = 3
    ),
    cex = 1.2
  );

  if (!legend) {
    custom_key = NULL;
  }
  
  if (!is.null(mappa)) {
    xlab = expression(paste(Log[10]," Mappable Locus Length"))
  } else {
    xlab = expression(paste(Log[10]," Locus Length"))
  }
	
	xmin_nopad = base::ifelse(is.null(xlim[1]),floor(min(gpw$log10_length)),floor(xlim[1]));
	xmax_nopad = base::ifelse(is.null(xlim[2]),ceiling(max(gpw$log10_length)),ceiling(xlim[2]));
	
	scales = list(
		x = list(
			axs = 'i',
			at = seq(xmin_nopad,xmax_nopad,1)
		),
		y = list(
			axs = 'i',
			at = seq(0,1,0.2)
		)
	);
	
  plotobj = xyplot(
    false_prob + prob_peak ~ log10_length,
    gpw,
    xlab=list(label=xlab,cex=1.4),
    ylab=list(label="P(peak)",cex=1.4),
    ylim=c(-0.05,1.05),
		xlim=c(xmin_nopad - 0.5,xmax_nopad + 0.5),
    panel=panel_func,
    type="l",
    lwd=3,
    key=custom_key,
    scales=list(cex=1.4),
    par.settings=simpleTheme(pch=c(15,17),col=c(col_rand_peak,col_spline))
  ) + as.layer(xyplot(peak ~ log_avg_length,avg_bins,pch=20,cex=0.4,col="black"));
    
  #trellis.focus("panel",1,1,highlight=F);
  #gpw_peak = subset(gpw,peak == 1);
  #gpw_nopeak = subset(gpw,peak == 0);
  #panel.rug(gpw_peak$log10_length,regular=F,col="black");
  #panel.rug(gpw_nopeak$log10_length,regular=T,col="black");
  #trellis.unfocus();

  return(plotobj);
}

plot_spline_mappa = function(locusdef,peak_genes,mappa) {
  # Calculate spline for mappability. 
  gpw = calc_spline_mappa(locusdef,peak_genes,mappa);
  
  # Add in empirical peak by binning. 
  gpw = add_emp_peak_mappa(gpw);
  
  # Average peak/lengths. 
  avg_bins = avg_binned_peak_mappa(gpw,bin_size=25);
  
  # Order by mappability. 
  gpw = gpw[order(gpw$mappa),];
  
  panel_func = function(x,y,...) { 
    panel.xyplot(x,y,...); 
    panel.abline(h=mean(y),col="blue",lwd=3);
  }

  custom_key = list( 
    text = list(
      c(
        "Mean(peak)",
        "Spline Fit - Mappability Only",
        "Proportion of Peaks in Bin"
      )
    ),
    lines = list(
      pch = c(20,17,20),
      type = c("l","l","p"),
      col = c("blue","green","black")
    )
  );

  plotobj = xyplot(
    prob_peak ~ mappa,
    gpw,
    xlab=list(label="Mappability",cex=1.25),
    ylab=list(label="P(peak)",cex=1.25),
    ylim=c(-0.05,1.05),
    panel=panel_func,
    type="l",
    lwd=3,
    key=custom_key,
    scales=list(cex=1.1),
    par.settings=simpleTheme(pch=c(17),col=c("green")),
  );

  trellis.focus("panel",1,1,highlight=F);
  grid.points(avg_bins$avg_mappa,avg_bins$peak,pch=20,gp=gpar(cex=0.35));
  trellis.unfocus();

   return(plotobj);
}

# For each peak, find the nearest TSS, and distance to it. 
peak_nearest_tss = function(peaks,tss,midpoint=T) {
  results = list();
  for (chrom in names(peaks)) {
    if (!chrom %in% names(tss$ranges)) {
      next;
    }
    
    peak_ranges = peaks[[chrom]];
    
    if (midpoint) { 
      peak_mids = mid(peak_ranges);
      peak_ranges = IRanges(peak_mids,peak_mids);
    }
    
    tss_ranges = tss$ranges[[chrom]];
    
    nearest_tss = tss_ranges[nearest(peak_ranges,tss_ranges)];
    
    dist_to_tss = distance(peak_ranges,nearest_tss);
    
		if (midpoint) {
			d = data.frame(
					chrom=chrom,
					peak_midpoint=start(peak_ranges),
					nearest_tss=start(nearest_tss),
					dist_to_tss=dist_to_tss,
					geneid=names(nearest_tss)
			);
		} else {
			d = data.frame(
					chrom=chrom,
					peak_start=start(peak_ranges),
					peak_end=end(peak_ranges),
					nearest_tss=start(nearest_tss),
					dist_to_tss=dist_to_tss,
					geneid=names(nearest_tss)
			);
		}
    
    results[[chrom]] = d;
  }
  
  return(merge_all(results));
}

plot_dist_to_tss = function(peaks,genome='hg19') {
	# Get peaks from user's file. 
	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks);
	} else if (class(peaks) == "character") {
    if (get_ext(peaks) == "bed") {
      message("Reading BED file: ",peaks);
      peakobj = read_bed(peaks);
    } else {
      message("Reading peaks file: ",peaks);
      peakobj = read_peaks(peaks);
    }
	}
	
	# Load TSS site info. 
  tss_code = sprintf("tss.%s",genome);
  data(list=tss_code,package = "chipenrich.data");
  tss = get(tss_code);
	
	plotobj = ..plot_dist_to_tss(peakobj,tss);
	return(plotobj);
}

..plot_dist_to_tss = function(peaks,tss) {
  # Calculate distance to each TSS. 
  tss_peaks = peak_nearest_tss(peaks,tss);
  
  # Create distance to TSS plot.
  max_dist = max(tss_peaks$dist_to_tss);
  breaks = breaks=c(0,100,1000,5000,10000,50000,100000,max_dist);
  dist_table = table(cut(tss_peaks$dist_to_tss,breaks=breaks));
  dist_table = dist_table / sum(dist_table);
  
  names(dist_table) = c("< 0.1","0.1 - 1","1 - 5","5 - 10","10 - 50","50 - 100","> 100")
  
  pf = function(...) { 
    args <- list(...);
    bar_labels = sprintf("%0.1f%%",args$y * 100);
    panel.text(seq(1,length(args$x)),args$y + 0.03,bar_labels,cex=1.5);
    panel.barchart(...);
  }
  
  plotobj = barchart(
    dist_table,
    panel=pf,
    horizontal=F,
    scales=list(rot=45,cex=1.6),
    col="gray",
    ylim=c(0,1),
    ylab=list(label="% of Peaks",cex=1.65),
    xlab=list(label="Distance to TSS (kb)",cex=1.65),
    main=list(label="Distribution of Distance from Peaks to Nearest TSS",cex=1.45)
  );
  
  return(plotobj);
}

genome_to_organism = function(genome) {
  code = substr(as.character(genome),1,2);
  if (code == 'mm') {
    org = 'mmu';
  } else if (code == 'hg') {
    org = 'hsa';
  } else if (code == 'rn') {
    org = 'rno';
  }
  else {
    org = NULL;
  }
  
  if (is.null(org)) {
    stop("Error: genome requested is not supported.");
  }
  
  org;
}

supported_locusdefs = function() {
  piqr = data(package = "chipenrich.data");
  data_files = piqr$results[,3];
  
  ldefs = grep("locusdef",data_files,value=T);
  unique(str_replace(ldefs,"locusdef\\.(.+?)\\.",""));
}

supported_read_lengths = function() {
  piqr = data(package = "chipenrich.data");
  data_files = piqr$results[,3];
  
  mappas = grep("mappa",data_files,value=T);
  sort(unique(as.numeric(str_replace(str_extract(mappas,"(\\d+)mer"),"mer",""))));
}

supported_genesets = function() {
  piqr = data(package = "chipenrich.data");
  data_files = piqr$results[,3];
  
  geneset_files = grep("geneset",data_files,value=T);
  unlist(unique(Map(function(x) x[[2]],str_split(geneset_files,"\\."))));
}

supported_genomes = function() {
  piqr = data(package = "chipenrich.data");
  data_files = piqr$results[,3];

  na.omit(unique(str_match(data_files,"(locusdef|mappa)\\.(\\w+)\\.")[,3]));
}

supported_methods = function() {
  return(names(SUPPORTED_METHODS));
}

# Checks to see if all elements of a list argument are possible. 
check_arg = function(arg,possible_values,value=F) {
  if (!all(arg %in% possible_values)) {
    if (value) {
      return(arg[!arg %in% possible_values]);
    } else {
      return(F);
    }
  } else {
    return(T);
  }
}

read_mappa = function(file_path) {
  if (!file.exists(file_path)) {
    stop("Error: could not find mappability file: ",file_path);
  }

	d = read.table(file_path,sep="\t",header=T,stringsAsFactors=F);
	
	# Check columns. 
	for (col in c("geneid","mappa")) {
		if (!col %in% names(d)) {
			stop(sprintf("Error reading mappability data: no column named '%s' in file.",col));
		}
	}
	
	# Genes in this file should not be duplicated. 
	if (sum(duplicated(d$geneid)) > 0) {
		stop("Error reading mappability data: duplicate geneids exist in file.");
	}
  
  # Mappability should be between 0 and 1. 
	if (min(d$mappa) < 0) {
    stop("Error: mappability must be >= 0.");
  }
  
  if (max(d$mappa) > 1) {
    stop("Error: mappability must be <= 1.");
  }
  
	return(d);
}

filter_genesets = function(x,max_geneset_size = 2000) {
  x_class = class(x);
  x_class_attr = attr(x_class,"package");
  
  if (is.null(x_class) || is.null(x_class_attr)) {
    stop("Error: bad geneset object in filtering step..");
  }
  
  if (x_class != "GeneSet" || x_class_attr != "chipenrich.data") {
    stop("Error: bad geneset object in filtering step..");
  }
  
  g = as.list(x@set.gene);
  g = Filter(function(x) length(x) <= max_geneset_size,g);
  x@set.gene = as.environment(g);
  
  return(x);
}

# Recodes # of peaks to be: 
# { 0, if no peaks
# { 1, if number of peaks is >= threshold argument
recode_peaks = function(num_peaks,threshold=1) {
  as.numeric(num_peaks >= threshold);
}

# Creates a LocusDefinition object given a flat file of definitions. 
setup_ldef = function(filepath) {
  # Check if file exists. 
  if (!file.exists(filepath)) {
    stop("Error: unable to open file for reading: ",filepath);
  } 

  # Read in flat file. 
  message("Reading user-specified gene locus definitions: ",filepath);
  d = read.table(filepath,sep="\t",header=T,stringsAsFactors=F);
  
  # Warn if missing entries, we'll ignore them though. 
  num_missing = sum(!complete.cases(d));
  if (num_missing > 0) {
    message(sprintf("Warning: %i rows of user-provided locus definition were missing values, skipping these rows..",num_missing));
  }
  
  # Remove rows with missing data. 
  d = na.omit(d);
  
  filename = strip_ext(filepath);
  
  # Create new locus definition object. 
  object = new("LocusDefinition");

  # Remove duplicated rows, not duplicated genes.
  # Some genes will exist on multiple rows because they are split by
  # other transcripts, or small nuclear RNAs, etc.
  d = unique(d);

  # Create an IRanges object representing the loci for each gene on that chromosome. 
  chroms = list();
  for (chr in unique(d$chrom)) {
    genes_chrom = subset(d,chrom == chr);
    chroms[[chr]] = IRanges(start=genes_chrom$start,end=genes_chrom$end,names=genes_chrom$geneid);
  }

  object@dframe = d;
  object@chrom2iranges = chroms;

  # Store as GRanges object as well, for convenience.
  object@granges = GRanges(
    seqnames=d$chrom,
    ranges=IRanges(d$start,d$end),
    names=d$geneid
  );
  
  # Check to make sure locus definitions are disjoint - that is, they do not overlap each other. 
  if(!isDisjoint(object@granges)) {
    stop("Error: user-provided locus definitions overlap - there should be disjoint ranges for all genes.");
  }

  object;
}

# d_fwf = function(d) {
  # d = colwise(as.character)(d);
  # widths = sapply(d,function(x) max(nchar(x)));
  # widths + 2;
  
  # final = NULL;
  # for (i in 1:dim(d)[2]) {
    # n = names(d)[i];
    # sf = sprintf("%%-0%is",widths[i]);
    # v = sprintf(sf,d[,n]);
    # final = cbind(final,v);
  # }
  
  # final = data.frame(final,stringsAsFactors=F);
  # names(final) = names(d);
  
  # final;
# }

chipenrich = function(
  peaks,
  out_name = "chipenrich",
  out_path = getwd(),
  genome = "hg19",
  genesets = c(
    'GOBP',
    'GOCC',
    'GOMF'
  ),
  locusdef = "nearest_tss",
  method = "chipenrich",
  fisher_alt = "two.sided",
  use_mappability = F,
  mappa_file = NULL,
  read_length = 36,
  qc_plots = T,
  max_geneset_size = 2000,
  num_peak_threshold = 1
) {

  l = unlist(as.list(environment()));
  opts = data.frame(
    args = names(l),
    values = sapply(unlist(l),paste,collapse=","),
    stringsAsFactors = F
  );
  rownames(opts) = 1:length(l);
  
  # Check genome. 
  if (!genome %in% supported_genomes()) {
    stop("genome not supported: ",genome);
  }
  
  # Organism - human or mouse. 
  organism = genome_to_organism(genome);
  
  # Check genesets. The user can provide multiple genesets, and 
  # each one will be tested. 
  if (!check_arg(genesets,supported_genesets())) {
    bad_args = check_arg(genesets,supported_genesets(),value=T);
    stop("bad geneset(s) requested: ",paste(bad_args,collapse=", "));
  }
  
  # Check and load locus definition. 
  user_defined_ldef = file.exists(locusdef);
  if (user_defined_ldef) {
    # Load user-defined locus definition file. 
    ldef = setup_ldef(locusdef);
  } else {
    if (!locusdef %in% supported_locusdefs()) {
      stop("Error: invalid definition requested: ",locusdef);
    }
    
    # Load locus definitions. 
    ldef_code = sprintf("locusdef.%s.%s",genome,locusdef);
    data(list=ldef_code,package = "chipenrich.data");
    ldef = get(ldef_code);
  }
  
  # If the user specified their own mappability, we can't use the built in mappability - they must provide their own. 
  user_defined_mappa = !is.null(mappa_file);
  
  # The user specified mappability - if they accidentally also thought use_mappability had to be enabled, just disable it. 
  if (user_defined_mappa & use_mappability) {
    use_mappability = FALSE;
  }
  
  if (user_defined_ldef) {
    if (use_mappability) {
      message("Warning: built-in mappability cannot be used with a user-defined locus definition, you must calculate mappability for your definition and pass it in with the mappa_file argument.");
      use_mappability = FALSE;
    }
  }
  
  # Check read length. 
  if (use_mappability) {
    if (!as.numeric(read_length) %in% supported_read_lengths()) {
      stop("Error: bad read length requested: ",read_length);
    }
  }
  
  # Load mappability if requested. 
  if (user_defined_mappa) {
    message("Reading user-specified gene locus mappability file: ",mappa_file);
    mappa = read_mappa(mappa_file);
  } else {
    if (use_mappability) {
      mappa_code = sprintf("mappa.%s.%s.%imer",genome,locusdef,read_length);
      data(list=mappa_code,package = "chipenrich.data");
      mappa = get(mappa_code);
      mappa = na.omit(mappa);
    } else {
      mappa = NULL;
    }
  }
  
  # If they specified their own mappability, we should check that the gene names overlap. 
  if (user_defined_mappa && user_defined_ldef) {
    total_unique_genes = union(mappa$geneid,ldef@dframe$geneid);
    mappa_ldef_inters = intersect(mappa$geneid,ldef@dframe$geneid);
    frac_overlap = length(mappa_ldef_inters) / length(total_unique_genes);
    
    if (frac_overlap < 0.95) { 
      stop("Error: your mappability genes and locus definition genes overlap by less than 95% (they should match almost exactly)..");
    }
  }
  
  # Testing function. 
  get_test_method = function(x) {
    if (method %in% names(SUPPORTED_METHODS)) {
      return(SUPPORTED_METHODS[[method]]);
    } else if (method %in% names(HIDDEN_METHODS)) {
      return(HIDDEN_METHODS[[method]]);
    } else {
      stop(sprintf("Error: invalid enrichment test requested: %s, contact developer.",method));
    }
  }
  testf = get_test_method(method);
  test_func = get(testf);
  
	# Test name. 
	method_name = METHOD_NAMES[[method]];
 
  # Get peaks from user's file. 
	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks);
	} else if (class(peaks) == "character") {
    if (str_sub(peaks,-7,-1) == ".bed.gz") {
      message("Reading BED file: ",peaks);
      peakobj = read_bed(peaks);
    } else if (str_sub(peaks,-4,-1) == ".bed") {
      message("Reading BED file: ",peaks);
      peakobj = read_bed(peaks);
    } else {
      message("Reading peaks file: ",peaks);
      peakobj = read_peaks(peaks);
    }
	}
  
  # Warn user if they are trying to use FET with a 
  # locus definition that might lead to biased results. 
  #if ((method == "fet") & (!locusdef %in% c("1kb","5kb"))) {
  if (method == "fet") {
    if (is.character(locusdef) && !locusdef %in% c("1kb","5kb")) {
      message("");
      message("Warning: Fisher's exact test should only be used with the 1kb or 5kb locus definition.");
      message("");
    }
    
    if (user_defined_ldef) {
      message("");
      message("Warning: Fisher's exact test may give biased results if the spline fit for the gene locus definitions is not flat (see QC plots.)");
      message("");
    }
  }
  
  # Warn user if they are using the binomial test. 
  if (method == "binomial") {
    message("");
    message("Warning: the binomial test is provided for comparison purposes only.");
    message("This test will almost always give biased results favoring gene sets with short average locus length.");
    message("");
  }
  
	# Number of peaks in data. 
  num_peaks = sum(sapply(peakobj,function(x) length(x)))
  
  # Load genesets and filter them if necessary. 
  geneset_list = list();
  for (gs in genesets) {
    geneset_code = sprintf("geneset.%s.%s",gs,organism);
    data(list=geneset_code,package = "chipenrich.data");
    
    geneset_list[[geneset_code]] = filter_genesets(get(geneset_code),max_geneset_size);
  }
  
  # Load TSS site info. 
  tss_code = sprintf("tss.%s",genome);
  data(list=tss_code,package = "chipenrich.data");
  tss = get(tss_code);
  
  # Assign peaks to genes. 
  message("Assigning peaks to genes..");
  assigned_peaks = assign_peaks(peakobj,ldef,tss,midpoint=T);
  peak_genes = unique(assigned_peaks$geneid);
  	
	# Add gene symbols to peak genes. 
	genes_code = sprintf("genes.%s",organism);
	data(list=genes_code,package = "chipenrich.data");
	gene2symbol = get(genes_code);
	gene2symbol = change_names(gene2symbol,list(GENEID="geneid",SYMBOL="gene_symbol"));
	assigned_peaks = merge(assigned_peaks,gene2symbol,by="geneid",all.x=T);
  #  geneid chrom peak_start  peak_end peak_midpoint gene_locus_start gene_locus_end nearest_tss nearest_tss_gene dist_to_tss nearest_tss_gene_strand gene_symbol
	column_order = c(
		"chrom",
		"peak_start",
		"peak_end",
		"peak_midpoint",
		"geneid",
		"gene_symbol",
    "gene_locus_start",
    "gene_locus_end",
    "nearest_tss",
    "dist_to_tss",
    "nearest_tss_gene",
    "nearest_tss_gene_strand"
	);
	assigned_peaks = assigned_peaks[,column_order];
  
  ppg = num_peaks_per_gene(assigned_peaks,ldef,mappa);
  ppg$peak = recode_peaks(ppg$num_peaks,num_peak_threshold);
  	  
  # Run chipenrich method on each geneset. 
  results = list();
  for (gobj in geneset_list) {
		message(sprintf("Test: %s",method_name));
		message(sprintf("Genesets: %s",gobj@type));
		message("Running tests..");
    if (testf == "test_gam") {
      rtemp = test_func(gobj,ppg);
    }
    if (testf == "test_fisher_exact") {
      rtemp = test_func(gobj,ppg,alternative=fisher_alt);
    }
    if (testf == "test_binomial") {
      rtemp = test_func(gobj,ppg);
    }
    
    # Annotate with geneset descriptions. 
    rtemp$"Description" = as.character(mget(rtemp$Geneset.ID,gobj@set.name,ifnotfound=NA));
    rtemp$"Geneset.Type" = gobj@type;
    
    results[[gobj@type]] = rtemp;
  }
  enrich = merge_all(results);
  
  # Re-order the columns to something sensible. 
	column_order = c(
    "Geneset.Type",
    "Geneset.ID",
    "Description",
    "P.value",
    "FDR",
    "Effect",
		"Odds.Ratio",
		"P.Success",
    "Status",
    "N.Geneset.Genes",
    "N.Geneset.Peak.Genes",
		"Geneset.Peak.Genes"
  );
	column_order = intersect(column_order,names(enrich));
  enrich = enrich[,column_order];
  
  # Order results by p-value. 
  enrich = enrich[order(enrich$P.value),];
  
  # If there is a status column, re-sort so enriched terms are on top. 
  if ("Status" %in% names(enrich)) {
    enrich = enrich[order(enrich$Status,decreasing=T),];
  }
  
  # Pull out tests that failed. 
  bad_enrich = subset(enrich,is.na(P.value));
  enrich = subset(enrich,!is.na(P.value));

  # Write results to file. 
  if (!is.null(out_name)) {
    filename_analysis = file.path(out_path,sprintf("%s_results.tab",out_name));
    write.table(enrich,file=filename_analysis,row.names=F,quote=F,sep="\t");
    message("Wrote results to: ",filename_analysis);
		
    filename_peaks = file.path(out_path,sprintf("%s_peaks.tab",out_name));
		write.table(assigned_peaks,file=filename_peaks,row.names=F,quote=F,sep="\t");
		message("Wrote peak-to-gene assignments to: ",filename_peaks);
    
    filename_opts = file.path(out_path,sprintf("%s_opts.tab",out_name));
		write.table(opts,file=filename_opts,row.names=F,quote=F,sep="\t");
		message("Wrote run options/arguments to: ",filename_opts);
    
    filename_ppg = file.path(out_path,sprintf("%s_peaks-per-gene.tab",out_name));
		write.table(ppg,file=filename_ppg,row.names=F,quote=F,sep="\t");
		message("Wrote count of peaks per gene to: ",filename_ppg);
    
    # If the user requested QC plots, generate those as well. 
    # -- Spline fit plot
    # -- Histogram of distance from peaks to TSSs
    # -- Histogram of p-values from test
    # -- Expected # peaks vs. observed # peaks
    if (qc_plots) {
      filename_qcplots = file.path(out_path,sprintf("%s_qcplots.pdf",out_name));
      pdf(filename_qcplots);
      print(..plot_spline_length(ldef,peak_genes,num_peaks,mappa=mappa));
      print(..plot_dist_to_tss(peakobj,tss));
#      print(..plot_expected_peaks(ppg));
      dev.off();
      
      message("Wrote QC plots to: ",filename_qcplots);
    }
  }

  return(list(
    peaks = assigned_peaks,
    results = enrich,
    opts = opts,
    peaks_per_gene = ppg
  ));
}
