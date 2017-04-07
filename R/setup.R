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

#' Function to read custom mappability files
#'
#' This function reads a two-columned tab-delimited text file (with header). Expected
#' column names are 'mappa' and 'gene_id'. Each line is for a unique 'gene_id'
#' and contains the mappability (between 0 and 1) for that gene.
#'
#' @param file_path A file path for the custom mappability.
#'
#' @return A \code{data.frame} containing \code{gene_id} and \code{mappa} columns.
read_mappa = function(file_path) {

	d = read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

	# Check columns.
	if(!(all(colnames(d) %in% c('mappa','gene_id')))) {
		stop("Error: header must contain columns named 'mappa' and 'gene_id'")
	}

	# Genes in this file should not be duplicated.
	if (sum(duplicated(d$gene_id)) > 0) {
		stop("Error reading mappability data: duplicate gene_ids exist in file.")
	}

	# Mappability should be between 0 and 1.
	if (min(d$mappa) < 0) {
		stop("Error: mappability must be >= 0.")
	}

	if (max(d$mappa) > 1) {
		stop("Error: mappability must be <= 1.")
	}

	return(d)
}

#' Function to read custom locus definition from file
#'
#' This function reads a tab-delimited text (with a header) file that should have
#' columns 'chr', 'start', 'end', and a column named 'gene_id' (or 'geneid') with the
#' Entrez Gene ID. If a \code{supported_genomes()} is given, then a column
#' of gene symbols named 'symbol', will be added. If an unsupported genome is
#' used there are two options: 1) Have a column named 'symbols' with the gene
#' symbols in the custom locus definition, and leave \code{genome = NA}, or 2)
#' leave \code{genome = NA}, do not provide gene symbols, and NAs will be used.
#'
#' @param file_path A file path for the custom locus definition.
#' @param genome A genome from \code{supported_genomes()}, default \code{NA}.
#'
#' @return A \code{LocusDefinition} class object with slots \code{dframe},
#' \code{granges}, \code{genome.build}, and \code{organism}.
setup_ldef = function(file_path, genome = NA) {

	message("Reading user-specified gene locus definitions: ", file_path)
	df = read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

	# Remove duplicated rows, not duplicated genes.
	# Some genes will exist on multiple rows because they are split by
	# other transcripts, or small nuclear RNAs, etc.
	df = unique(df)

	# Rename 'geneid' columns 'gene_id' for backwards compatibility
	if('geneid' %in% colnames(df)) {
		colnames(df)[which(colnames(df) == 'geneid')] = 'gene_id'
	}

	if(!('gene_id' %in% colnames(df))) {
		stop("Error: Custom locus definition must have column named 'gene_id'.")
	}

	if('symbol' %in% colnames(df)) {
		# Do nothing
		message('Using given symbol column...')
	} else if (genome %in% supported_genomes()) {
		message('Using orgDb package to get gene symbols...')
		eg2symbol = genome_to_orgdb(genome)
		df$symbol = eg2symbol[match(df$gene_id, eg2symbol$gene_id), 'symbol']
	} else {
		message('Setting gene symbols to NA...')
		df$symbol = NA
	}

	if(genome %in% supported_genomes()) {
		gr = GenomicRanges::makeGRangesFromDataFrame(
			df = df,
			seqinfo = GenomeInfoDb::Seqinfo(genome = genome),
			keep.extra.columns = TRUE)
	} else {
		gr = GenomicRanges::makeGRangesFromDataFrame(
			df = df,
			keep.extra.columns = TRUE)
	}

	# Create new locus definition object.
	object = new("LocusDefinition")

	object@dframe = df
	object@granges = gr

	return(object)
}

# Creates an object that mimics the GeneSet class
# from the chipenrich.data package, as in:
#    setClass("GeneSet",representation(
#      set.gene = "environment",
#      type = "character",
#      set.name = "environment",
#      all.genes = "character",
#      organism = "character",
#      dburl = "character"
#    ),prototype(
#      set.gene = new.env(parent=emptyenv()),
#      type = "",
#      set.name = new.env(parent=emptyenv()),
#      all.genes = "",
#      organism = "",
#      dburl = ""
#    ),
#      package = "chipenrich.data"
#    )

#' Function read custom gene sets from file
#'
#' This function reads a two-columned tab-delimited text file (with header). Column
#' names are ignored, but the first column should be geneset names or IDs and the
#' second column should be Entrez Gene IDs.
#'
#' @param file_path A file path for the custom gene set.
#'
#' @return A \code{GeneSet} class object.
setup_geneset = function(file_path) {

	# Read in flat file.
	message("Reading user-specified gene set definitions: ", file_path)
	d = read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

	# Split the gene_id column (second) by the geneset_id column (first)
	gs = split(d[,2], d[,1])

	# Create the names list which is identity
	gs.names = as.list(names(gs))
	names(gs.names) = names(gs)

	# Create new GeneSet object
	object = new('GeneSet')

	# Populate the GeneSet object
	object@type = 'user-supplied'
	object@organism = 'user-supplied'
	object@dburl = 'user-supplied'

	object@set.gene = as.environment(gs)
	object@all.genes = as.character(unlist(gs, use.names = FALSE))

	object@set.name = as.environment(gs.names)
	message('Done setting up user-specified geneset..')

	return(object)
}
