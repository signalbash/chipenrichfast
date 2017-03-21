#' chipenrich: Gene Set Enrichment For ChIP-seq Peak Data
#'
#' ChIP-Enrich performs gene set enrichment testing using peaks called from a
#' ChIP-seq experiment. The method empirically corrects for confounding factors
#' such as the length of genes, and the mappability of the sequence surrounding genes.
#'
#' @docType package
#' @name chipenrich
#'
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
#' @import parallel
#' @importFrom plyr rbind.fill
#' @import rms
#' @import rtracklayer
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom stats as.formula binom.test coef complete.cases fisher.test fitted na.omit p.adjust pchisq predict qpois quantile resid
#' @importFrom stringr str_match str_sub
#' @importFrom utils data read.table tail write.table
NULL
