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

#' Function to setup locus definitions
#'
#' @param ldef_code One of 'nearest_tss', 'nearest_gene', 'exon', 'intron', '1kb',
#' '1kb_outside', '1kb_outside_upstream', '5kb', '5kb_outside', '5kb_outside_upstream',
#' '10kb', '10kb_outside', '10kb_outside_upstream'. Alternately, a file path for
#' a custom locus definition. NOTE: Must be for a \code{supported_genome()}, and
#' must have columns 'chr', 'start', 'end', and 'gene_id', or 'geneid'.
#' @param genome One of the \code{supported_genomes()}.
#' @param randomization One of \code{NULL}, 'complete', 'bylength', or 'bylocation'.
#' See the Randomizations section in \code{?chipenrich}. Default NULL.
#'
#' @return A list with components \code{ldef} and \code{tss}.
#' @include randomize.R
setup_locusdef = function(ldef_code, genome, randomization = NULL) {
	user_defined_ldef = file.exists(ldef_code)

	if (user_defined_ldef) {
		ldef = read_ldef(ldef_code, genome)
	} else {
		if (!any(supported_locusdefs()$genome == genome & ldef_code == supported_locusdefs()$locusdef)) {
			stop(sprintf("Error: invalid genome / definition combination requested: %s %s", genome, ldef_code))
		}

		ldef_rdata_code = sprintf("locusdef.%s.%s", genome, ldef_code)
		data(list = ldef_rdata_code, package = "chipenrich.data", envir = environment())
		ldef = get(ldef_rdata_code)
	}

	# Load TSS site info.
	tss_code = sprintf("tss.%s", genome)
	data(list=tss_code, package = "chipenrich.data", envir = environment())
	tss = get(tss_code)

	# Deal with randomizations
	if(is.null(randomization)) {

	} else if (randomization == 'complete') {
		ldef = randomize_ldef_complete(ldef)
	} else if (randomization == 'bylength') {
		ldef = randomize_ldef_bylength(ldef, resolution = 100)
	} else if (randomization == 'bylocation') {
		ldef = randomize_ldef_bylocation(ldef, resolution = 50)
	}

	return(list(
		ldef = ldef,
		tss = tss))
}

#' Function to setup genesets
#'
#' @param gs_codes A character vector of geneset databases to be tested for
#' enrichment. See \code{supported_genesets()}. Alternately, a file path to a
#' a tab-delimited text file with header and first column being the geneset ID
#' or name, and the second column being Entrez Gene IDs.
#' @param ldef_obj A \code{LocusDefinition} object to use for filtering gene sets
#' based on which genes are defined in the locus defintion.
#' @param genome One of the \code{supported_genomes()}.
#' @param min_geneset_size Sets the minimum number of genes a gene set may have
#' to be considered for enrichment testing.
#' @param max_geneset_size Sets the maximum number of genes a gene set may have
#' to be considered for enrichment testing.
#'
#' @return A list with components consisting of \code{GeneSet} objects for each
#' of the elements of \code{genesets}. NOTE: Custom genesets must be run separately
#' from built in gene sets.
setup_genesets = function(gs_codes, ldef_obj, genome, min_geneset_size, max_geneset_size) {
	if(class(ldef_obj) != 'LocusDefinition') {
		stop('ldef_obj not of class LocusDefinition.')
	}

	organism = genome_to_organism(genome)

	user_defined_geneset = file.exists(gs_codes)

	if(user_defined_geneset) {
		geneset_list = list()
		geneset_code = 'user-supplied'

		geneset_list[[geneset_code]] = read_geneset(gs_codes)
		geneset_list[[geneset_code]] = filter_genesets(geneset_list[[geneset_code]], ldef_obj, min_geneset_size, max_geneset_size)
	} else {
		if (!any(supported_genesets()$organism == organism & gs_codes %in% supported_genesets()$geneset)) {
			stop("Invalid organism / geneset combination requested. Please see supported_genesets().")
		}

		# If the user does not provide a path to the geneset.
		# That is, do the normal thing
		geneset_list = list()
		for (gs_code in gs_codes) {
			gs_rdata_code = sprintf("geneset.%s.%s", gs_code, organism)
			data(list = gs_rdata_code, package = "chipenrich.data", envir = environment())

			geneset_list[[gs_rdata_code]] = filter_genesets(get(gs_rdata_code), ldef_obj, min_geneset_size, max_geneset_size)
		}
	}

	return(geneset_list)
}

#' Function to setup mappability
#'
#' @param mappa_code One of \code{NULL}, a file path to a custom mappability file,
#' or an \code{integer} for a valid read length given by \code{supported_read_lengths}.
#' If a file, it should contain a header with two column named 'gene_id' and 'mappa'.
#' Gene IDs should be Entrez IDs, and mappability values should range from 0 and 1.
#' Default value is NULL.
#' @param genome One of the \code{supported_genomes()}.
#' @param ldef_code One of 'nearest_tss', 'nearest_gene', 'exon', 'intron', '1kb',
#' '1kb_outside', '1kb_outside_upstream', '5kb', '5kb_outside', '5kb_outside_upstream',
#' '10kb', '10kb_outside', '10kb_outside_upstream'. Alternately, a file path for
#' a custom locus definition. NOTE: Must be for a \code{supported_genome()}, and
#' must have columns 'chr', 'start', 'end', and 'gene_id', or 'geneid'.
#' @param ldef_obj A \code{LocusDefinition} object.
#'
#' @return A \code{data.frame} with columns \code{gene_id} and \code{mappa}.
setup_mappa = function(mappa_code, genome, ldef_code, ldef_obj) {
	if(is.null(mappa_code)) {
		return(NULL)
	}

	if(class(ldef_obj) != 'LocusDefinition') {
		stop('ldef_obj not of class LocusDefinition.')
	}

	# mappa_code		ldef_code		action
	# NULL				file			return(NULL)
	# NULL				built-in		return(NULL)
	# file				file			check overlap
	# file				built-in		check overlap
	# read_length		file			stop
	# read_length		built-in		check valid mappa_code

	user_defined_mappa = file.exists(mappa_code)
	user_defined_ldef = file.exists(ldef_code)

	if(user_defined_mappa) {
		message("Reading user-specified gene locus mappability file: ", mappa_code)
		mappa = read_mappa(mappa_code)

		# Check overlap of user defined mappability gene_id with locus definition gene_id
		total_unique_genes = union(mappa$gene_id, ldef_obj@dframe$gene_id)
		mappa_ldef_inters = intersect(mappa$gene_id, ldef_obj@dframe$gene_id)
		frac_overlap = length(mappa_ldef_inters) / length(total_unique_genes)

		if (frac_overlap < 0.95) {
			stop("Error: your mappability genes and locus definition genes overlap by less than 95% (they should match almost exactly)..")
		}
	} else {
		if(user_defined_ldef) {
			stop("Built-in mappability cannot be used with a user-defined locus definition. Please calculate mappability for your definition and pass its file path to the 'mappability' parameter.")
		} else {
			if (!any(supported_read_lengths()$genome == genome & supported_read_lengths()$locusdef == ldef_code & supported_read_lengths()$read_length == as.numeric(mappa_code))) {
				stop(sprintf("Error: bad genome / locusdef / read length combination requested: %s %s %s", genome, ldef_code, mappa_code))
			}

			mappa_rdata_code = sprintf("mappa.%s.%s.%smer", genome, ldef_code, mappa_code)
			data(list = mappa_rdata_code, package = "chipenrich.data", envir = environment())
			mappa = get(mappa_rdata_code)
			mappa = na.omit(mappa)
			colnames(mappa) = c('gene_id', 'mappa')
		}
	}

	return(mappa)
}
