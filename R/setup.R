#' Function to filter genesets by locus definition and size
#'
#' This function filters gene sets based on the genes that are present in a particular
#' locus definition. After determining which genes are present in both the GeneSet,
#' \code{gs_obj}, and the LocusDefinition \code{ldef_obj}, gene sets are filtered
#' by size with \code{min_geneset_size} and \code{max_geneset_size}.
#'
#' @param gs_obj A valid GeneSet object
#' @param ldef_obj A valid LocusDefinition object
#' @param min_geneset_size An integer indicating the floor for genes in a geneset. Default 15.
#' @param max_geneset_size An integer indicating the ceiling for genes in a geneset. Default 2000.
#'
#' @return An altered \code{gs_obj} with changed \code{set.gene} and \code{all.genes}
#' slots reflecting \code{min_geneset_size} and \code{max_geneset_size} after intersecting
#' with the genes present in the particular locus definition.
filter_genesets = function(gs_obj, ldef_obj, min_geneset_size = 15, max_geneset_size = 2000) {
	if (class(gs_obj) != "GeneSet") {
		stop("Error: gs_obj not of class GeneSet")
	}

	if(class(ldef_obj) != 'LocusDefinition') {
		stop("Error: ldef_obj not of class LocusDefinition")
	}

	# Coerce to character, because gene sets store genes as character(), and
	# intersect() is much faster when the vectors are of the same type.
	ldef_gene_ids = as.character(unique(ldef_obj@granges$gene_id))
	tmp_set.gene = as.list(gs_obj@set.gene)

	tmp_set.gene = lapply(tmp_set.gene, function(gs){intersect(ldef_gene_ids, gs)})
	tmp_set.gene = Filter(function(gs) length(gs) >= min_geneset_size, tmp_set.gene)
	tmp_set.gene = Filter(function(gs) length(gs) <= max_geneset_size, tmp_set.gene)

	gs_obj@set.gene = as.environment(tmp_set.gene)
	gs_obj@all.genes = unique(as.character(unlist(tmp_set.gene, use.names = FALSE)))

	return(gs_obj)
}
