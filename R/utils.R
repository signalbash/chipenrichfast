# This is really a function to check for a valid method

# Recodes # of peaks to be:
# { 0, if no peaks
# { 1, if number of peaks is >= threshold argument
recode_peaks = function(num_peaks, threshold = 1) {
  as.numeric(num_peaks >= threshold);
}

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
