filter_genesets = function(x,max_geneset_size = 2000) {
	x_class = class(x)
	x_class_attr = attr(x_class, "package")

	if (is.null(x_class) || is.null(x_class_attr)) {
		stop("Error: bad geneset object in filtering step..")
	}

	if (x_class != "GeneSet" || x_class_attr != "chipenrich.data") {
		stop("Error: bad geneset object in filtering step..")
	}

	g = as.list(x@set.gene)
	g = Filter(function(x) length(x) <= max_geneset_size, g)
	x@set.gene = as.environment(g)

	return(x)
}

read_mappa = function(file_path) {
	if (!file.exists(file_path)) {
		stop("Error: could not find mappability file: ", file_path)
	}

	d = read.table(file_path, sep = "\t", header = T, stringsAsFactors = F)

	# Check columns.
	for (col in c("gene_id", "mappa")) {
		if (!col %in% names(d)) {
			stop(sprintf("Error reading mappability data: no column named '%s' in file.",col))
		}
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
#' @param filepath A valid file path for the custom locus definition.
#' @param genome A genome from \code{supported_genomes()}, default \code{NA}.
#'
#' @return A \code{LocusDefinition} class object with slots \code{dframe},
#' \code{granges}, \code{genome.build}, and \code{organism}.
setup_ldef = function(filepath, genome = NA) {
	if (!file.exists(filepath)) {
		stop("Error: unable to open file for reading: ", filepath)
	}

	message("Reading user-specified gene locus definitions: ", filepath)
	df = read.table(filepath, sep = "\t", header = T, stringsAsFactors = F)

	num_missing = sum(!complete.cases(df))
	if (num_missing > 0) {
		message(sprintf("Warning: %i rows of user-provided locus definition were missing values, skipping these rows.", num_missing))
	}
	df = na.omit(df)
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
# User supplied genesets are expected to come as two columned,
# tab-delimited text files with the first column being the geneset
# name and the second column being Entrez Gene ID
setup_geneset = function(filepath) {

	# Check if file exists.
	if (!file.exists(filepath)) {
		stop("Error: unable to open file for reading: ", filepath)
	}

	# Read in flat file.
	message("Reading user-specified gene set definitions: ", filepath)
	d = read.table(filepath, sep = "\t", header = T, stringsAsFactors = F)

	# Warn if missing entries, we'll ignore them though.
	num_missing = sum(!stats::complete.cases(d))
	if (num_missing > 0) {
		message(sprintf("Warning: %i rows of user-provided locus definition were missing values, skipping these rows..", num_missing))
	}

	# Remove rows with missing data.
	d = stats::na.omit(d)

	filename = strip_ext(filepath)

	# Create shell for geneset list
	# Use numerical column accessors
	gs.names = unique(d[,1])
	gs = as.list(gs.names)
	names(gs) = gs.names

	gs.names.list = gs

	# Populate the shell
	# Use numerical column accessors
	gs = lapply(gs, function(g){
		return(unique(subset(d, d[,1] == g)[,2]))
	})

	# Create new GeneSet object
	object = new('GeneSet')

	# Populate the GeneSet object
	object@type = 'user-supplied'
	object@organism = 'user-supplied'
	object@dburl = 'user-supplied'

	object@set.gene = as.environment(gs)
	object@all.genes = as.character(sort(Reduce(function(x,y) union(x,y), object@set.gene)))

	object@set.name = as.environment(gs.names.list)
	message('Done setting up user-specified geneset..')

	return(object)
}
