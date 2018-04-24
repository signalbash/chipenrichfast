#' Display supported locus definitions
#'
#' The locus definitions are defined as below. For advice on selecting a locus
#' definition, see the 'Selecting A Locus Definition' section below.
#' \describe{
#'   \item{nearest_tss:}{The locus is the region spanning the midpoints between the TSSs of adjacent genes.}
#'   \item{nearest_gene:}{The locus is the region spanning the midpoints between the boundaries of each gene, where a gene is defined as the region between the furthest upstream TSS and furthest downstream TES for that gene. If two gene loci overlap each other, we take the midpoint of the overlap as the boundary between the two loci. When a gene locus is completely nested within another, we create a disjoint set of 3 intervals, where the outermost gene is separated into 2 intervals broken apart at the endpoints of the nested gene.}
#'   \item{1kb:}{The locus is the region within 1kb of any of the TSSs belonging to a gene. If TSSs from two adjacent genes are within 2 kb of each other, we use the midpoint between the two TSSs as the boundary for the locus for each gene.}
#'   \item{1kb_outside_upstream:}{The locus is the region more than 1kb upstream from a TSS to the midpoint between the adjacent TSS.}
#'   \item{1kb_outside:}{The locus is the region more than 1kb upstream or downstream from a TSS to the midpoint between the adjacent TSS.}
#'   \item{5kb:}{The locus is the region within 5kb of any of the TSSs belonging to a gene. If TSSs from two adjacent genes are within 10 kb of each other, we use the midpoint between the two TSSs as the boundary for the locus for each gene.}
#'   \item{5kb_outside_upstream:}{The locus is the region more than 5kb upstream from a TSS to the midpoint between the adjacent TSS.}
#'   \item{5kb_outside:}{The locus is the region more than 5kb upstream or downstream from a TSS to the midpoint between the adjacent TSS.}
#'   \item{10kb:}{The locus is the region within 10kb of any of the TSSs belonging to a gene. If TSSs from two adjacent genes are within 20 kb of each other, we use the midpoint between the two TSSs as the boundary for the locus for each gene.}
#'   \item{10kb_outside_upstream:}{The locus is the region more than 10kb upstream from a TSS to the midpoint between the adjacent TSS.}
#'   \item{10kb_outside:}{The locus is the region more than 10kb upstream or downstream from a TSS to the midpoint between the adjacent TSS.}
#'   \item{exon:}{Each gene has multiple loci corresponding to its exons. Overlaps between different genes are allowed.}
#'   \item{intron:}{Each gene has multiple loci corresponding to its introns. Overlaps between different genes are allowed.}
#' }
#'
#' @section Selecting A Locus Definition:
#' For a transcription factor ChIP-seq experiment, selecting a particular locus
#' definition for use in enrichment testing can have implications relating to how
#' the TF regulates genes. For example, selecting the '1kb' locus definition will
#' imply that the biological processes found enriched are a result of TF regulation
#' near the promoter. In contrast, selecting the '5kb_outside' locus definition
#' will imply that the biological processes found enriched are a result of TF
#' regulation distal from the promoter.
#'
#' Selecting a locus definition can also help reduce the noise in the enrichment
#' tests. For example, if a TF is known to primarily regulate genes by binding
#' around the promoter, then selecting the '1kb' locus definition can help to
#' reduce the noise from TSS-distal peaks in the enrichment testing.
#'
#' The \code{\link{plot_dist_to_tss}} QC plot displays
#' where genomic regions fall relative to TSSs genome-wide, and can help inform
#' the choice of locus definition. For example, if many peaks fall far from the
#' TSS, the 'nearest_tss' locus definition may be a good choice because it will
#' capture all input genomic regions, whereas the '1kb' locus definition may
#' not capture many of the input genomic regions and adversely affect the
#' enrichment testing.
#'
#' @return A \code{data.frame} with columns \code{genome, locusdef}.
#'
#' @examples
#'
#' supported_locusdefs()
#'
#' @export
supported_locusdefs = function() {
	piqr = data(package = "chipenrich.data")
	data_files = piqr$results[,3]

	ldefs = grep("locusdef", data_files, value = TRUE)
	combos = Reduce(rbind,sapply(ldefs, strsplit, '[.]'))[,c(2,3)]
	df = data.frame(
		'genome' = combos[,1],
		'locusdef' = combos[,2],
		stringsAsFactors = FALSE)
	return(df)
}

#' Display supported read lengths for mappability
#'
#' @return A \code{data.frame} with columns \code{genome, locusdef, read_length}.
#'
#' @examples
#'
#' supported_read_lengths()
#'
#' @export
supported_read_lengths = function() {
	piqr = data(package = "chipenrich.data")
	data_files = piqr$results[,3]

	mappas = grep("mappa", data_files, value = TRUE)
	combos = Reduce(rbind,sapply(mappas, strsplit, '[.]'))[,c(2,3,4)]
	df = data.frame(
		'genome' = combos[,1],
		'locusdef' = combos[,2],
		'read_length' = gsub('mer','',combos[,3]),
		stringsAsFactors = FALSE)
	return(df)
}

#' Display supported genesets for gene set enrichment.
#'
#' @return A \code{data.frame} with columns \code{geneset, organism}.
#'
#' @examples
#'
#' supported_genesets()
#'
#' @export
supported_genesets = function(organism=NULL) {
	piqr = data(package = "chipenrich.data")
	data_files = piqr$results[,3]

	geneset_files = grep("geneset", data_files, value = TRUE)
	combos = Reduce(rbind,sapply(geneset_files, strsplit, '[.]'))[,c(2,3)]
	df = data.frame(
		'geneset' = combos[,1],
		'organism' = combos[,2],
		stringsAsFactors = FALSE)
	df = df[order(df$organism), ]
    
    if (!is.null(organism)){
        df = df[df$organism==organism,]
    }
    
	return(df)
}

#' Display supported genomes.
#'
#' @return A vector indicating supported genomes.
#'
#' @examples
#'
#' supported_genomes()
#'
#' @export
supported_genomes = function() {
	piqr = data(package = "chipenrich.data")
	data_files = piqr$results[,3]

	as.character(stats::na.omit(unique(stringr::str_match(data_files,"(locusdef|mappa)\\.(\\w+)\\.")[,3])))
}

#' Display supported gene set enrichment methods.
#'
#' @return A vector indicating supported methods for gene set enrichment.
#'
#' @examples
#'
#' supported_methods()
#'
#' @export
supported_methods = function() {
	return(SUPPORTED_COMBOS)
}
