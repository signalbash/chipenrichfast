#' chipenrich: Gene Set Enrichment For ChIP-seq Peak Data
#'
#' The \code{chipenrich} package includes three classes of methods that adjust for potential confounders of gene set enrichment testing (locus length and mappability of the sequence reads). The first, \code{chipenrich}, is designed for use with transcription-factor (TF) based ChIP-seq experiments and other DNA sequencing experiments with narrow genomic regions. The second, \code{polyenrich}, is similarly designed for TF based ChIP-seq, but where the number of peaks present in gene loci may be important. The third, \code{broadenrich}, is designed for use with histone modification based ChIP-seq experiments and other DNA sequencing experiments with broad genomic regions.
#'
#' @docType package
#' @name chipenrich_package
#'
#' @importFrom AnnotationDbi mappedkeys
#' @importFrom BiocGenerics unlist
#' @import chipenrich.data
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom grDevices dev.off pdf
#' @import grid
#' @importFrom IRanges IRanges
#' @import lattice
#' @import latticeExtra
#' @import methods
#' @import mgcv
#' @import org.Dm.eg.db
#' @import org.Dr.eg.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#' @import parallel
#' @importFrom plyr rbind.fill
#' @import rms
#' @import rtracklayer
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom stats as.formula binom.test coef complete.cases fisher.test fitted na.omit p.adjust pchisq predict qpois quantile resid
#' @importFrom stringr str_match str_sub
#' @importFrom utils data read.table tail write.table
NULL
