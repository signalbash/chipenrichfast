#' Get the correct organism code based on genome
#'
#' Data from \code{chipenrich.data} uses three letter organism codes for the
#' \code{GeneSet} objects. This function ensures the correct objects are loaded.
#'
#' @param genome One of the \code{supported_genomes()}.
#'
#' @return A string for the three letter organism code. Convention is first letter
#' of the first word in the binomial name, and first two letters of the second word
#' in the binomial name. 'Homo sapiens' is then 'hsa', for example.
genome_to_organism = function(genome = supported_genomes()) {
    genome = match.arg(genome)

    if (grepl('^mm', genome)) {
        org = 'mmu'
    } else if (grepl('^hg', genome)) {
        org = 'hsa'
    } else if (grepl('^rn', genome)) {
        org = 'rno'
    } else if (grepl('^dm', genome)) {
        org = 'dme'
    } else if (grepl('^danRer', genome)) {
        org = 'dre'
    }

    return(org)
}

#' Get Entrez ID to gene symbol mappings for custom locus definitions
#'
#' If a custom locus definition is one of the \code{supported_genomes()}, then
#' the gene symbol column of the custom locus definition is populated using the
#' appropriate orgDb package.
#'
#' @param genome One of the \code{supported_genomes()}.
#'
#' @return A \code{data.frame} with \code{gene_id} and \code{symbol} columns.
genome_to_orgdb = function(genome = supported_genomes()) {
    genome = match.arg(genome)

    if(genome %in% c('hg19','hg38')) {
        egSYMBOL = org.Hs.eg.db::org.Hs.egSYMBOL
    } else if (genome %in% c('mm9','mm10')) {
        egSYMBOL = org.Mm.eg.db::org.Mm.egSYMBOL
    } else if (genome %in% c('rn4','rn5','rn6')) {
        egSYMBOL = org.Rn.eg.db::org.Rn.egSYMBOL
    } else if (genome %in% c('dm3','dm6')) {
        egSYMBOL = org.Dm.eg.db::org.Dm.egSYMBOL
    } else if (genome == 'danRer10') {
        egSYMBOL = org.Dr.eg.db::org.Dr.egSYMBOL
    }

    ### Build Entrez ID to gene symbol mapping
        mapped_genes = AnnotationDbi::mappedkeys(egSYMBOL)
        eg2symbol = as.data.frame(egSYMBOL[mapped_genes])
        eg2symbol$gene_id = as.integer(eg2symbol$gene_id)

    return(eg2symbol)
}

#' Get the test function name from the method name
#'
#' The \code{method} comes from what is used in \code{chipenrich()} or in
#' \code{polyenrich()}.
#'
#' @param method A character for the method used. One of the \code{supported_methods}
#' or one of the \code{HIDDEN_METHODS} in \code{constants.R}.
#'
#' @return A singleton named character vector with value of the test function and
#' name of the method.
get_test_method = function(method) {
    if (method %in% names(SUPPORTED_METHODS)) {
        return(SUPPORTED_METHODS[[method]])
    } else if (method %in% names(HIDDEN_METHODS)) {
        return(HIDDEN_METHODS[[method]])
    } else {
        stop(sprintf("Error: invalid enrichment test requested: %s",method))
    }
}

#' Post process the \code{data.frame} of enrichment results
#'
#' @param enrich A \code{data.frame} of the enrichment results from \code{broadenrich()},
#' \code{chipenrich()}, or \code{polyenrich()} created by \code{rbind}ing the
#' list of enrichment results for each of the \code{genesets}.
#'
#' @return A reformatted \code{data.frame} with columns in a specific order, filtered
#' of enrichment tests that failed, and ordered first by enrichment 'Status' (if
#' present) and then 'P.value'.
post_process_enrichments = function(enrich) {
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
        "Geneset.Avg.Gene.Length",
        "Geneset.Avg.Gene.Coverage",
        "Geneset.Peak.Genes")
    column_order = intersect(column_order, names(enrich))
    enrich = enrich[, column_order]

    # Order results by p-value.
    enrich = enrich[order(enrich$P.value), ]

    # If there is a status column, re-sort so enriched terms are on top.
    if ("Status" %in% names(enrich)) {
        enrich = enrich[order(enrich$Status, decreasing=TRUE), ]
    }

    # Pull out tests that failed.
    bad_enrich = subset(enrich, is.na(enrich$P.value))
    enrich = subset(enrich, !is.na(enrich$P.value))
    rownames(enrich) = c(1:nrow(enrich))

    return(enrich)
}

#' Recode a vector of number of peaks to binary based on threshold
#'
#' @param num_peaks An \code{integer} vector representing numbers of peaks per gene.
#' @param threshold An \code{integer} specifying the minimum number of peaks
#' required to code as 1.
#'
#' @return An binary vector where an entry is 1 if the corresponding entry of
#' \code{num_peaks} is >= \code{threshold} and is otherwise 0.
recode_peaks = function(num_peaks, threshold = 1) {
  as.numeric(num_peaks >= threshold);
}

#' Reset n_cores for Windows
#'
#' We use \code{parallel::mclapply} for multicore geneset enrichment testing, but
#' this function supports more than one core if the OS is not Windows. If the OS
#' is windows, the number of cores (\code{mc.cores}) must be set to 1.
#'
#' @param n_cores An \code{integer} passed to \code{broadenrich()},
#' \code{chipenrich()}, or \code{polyenrich()} indicating the number of cores
#' to use for enrichment testing.
#'
#' @return Either the original \code{n_cores} if the OS is not Windows, or 1 if
#' the OS is Windows.
reset_ncores_for_windows = function(n_cores) {
    os = Sys.info()[1]
	if(os == 'Windows') {
		message('Setting n_cores = 1 because Windows detected as OS.')
		n_cores = 1
	}

    return(n_cores)
}
