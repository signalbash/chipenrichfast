## ------------------------------------------------------------------------
library(chipenrich)

## ---- warning = FALSE, message = FALSE-----------------------------------
data(peaks_E2F4, package = 'chipenrich.data')
data(peaks_H3K4me3_GM12878, package = 'chipenrich.data')

head(peaks_E2F4)
head(peaks_H3K4me3_GM12878)

## ---- echo=FALSE---------------------------------------------------------
#peaks_E2F4 = subset(peaks_E2F4, peaks_E2F4$chrom == 'chr1')
#peaks_H3K4me3_GM12878 = subset(peaks_H3K4me3_GM12878, peaks_H3K4me3_GM12878$chrom == 'chr1')

## ------------------------------------------------------------------------
supported_genomes()

## ------------------------------------------------------------------------
# Take head because it's long
head(supported_locusdefs())

## ---- eval=FALSE---------------------------------------------------------
#  chr	start	end	geneid
#  chr1	839460	839610	148398
#  chr1	840040	840190	148398
#  chr1	840040	840190	57801
#  chr1	840800	840950	148398
#  chr1	841160	841310	148398

## ------------------------------------------------------------------------
# Take head because it's long
head(supported_genesets())

## ---- eval=FALSE---------------------------------------------------------
#  gs_id	gene_id
#  GO:0006631	30
#  GO:0006631	31
#  GO:0006631	32
#  GO:0006631	33
#  GO:0006631	34
#  GO:0006631	35
#  GO:0006631	36
#  GO:0006631	37
#  GO:0006631	51
#  GO:0006631	131
#  GO:0006631	183
#  GO:0006631	207
#  GO:0006631	208
#  GO:0006631	215
#  GO:0006631	225

## ------------------------------------------------------------------------
# Take head because it's long
head(supported_read_lengths())

## ---- eval=FALSE---------------------------------------------------------
#  mappa	gene_id
#  0.8	8487
#  0.1	84
#  0.6	91
#  1	1000

## ---- warning = FALSE, message = FALSE-----------------------------------
gs_path = system.file('extdata','vignette_genesets.txt', package='chipenrich')
results = broadenrich(peaks = peaks_H3K4me3_GM12878, genome = 'hg19', genesets = gs_path,
	locusdef = "nearest_tss", qc_plots = FALSE, out_name = NULL, n_cores=1)
results.be = results$results
print(results.be[1:5,1:5])

## ---- warning = FALSE, message = FALSE-----------------------------------
# Without mappability
gs_path = system.file('extdata','vignette_genesets.txt', package='chipenrich')
results = chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path,
	locusdef = "nearest_tss", qc_plots = FALSE, out_name = NULL, n_cores = 1)
results.ce = results$results
print(results.ce[1:5,1:5])

## ---- warning = FALSE, message = FALSE-----------------------------------
# With mappability
gs_path = system.file('extdata','vignette_genesets.txt', package='chipenrich')
results = chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path,
	locusdef = "nearest_tss", mappability='24', qc_plots = FALSE,
	out_name = NULL,n_cores=1)
results.cem = results$results
print(results.cem[1:5,1:5])

## ---- warning = FALSE, message = FALSE-----------------------------------
gs_path = system.file('extdata','vignette_genesets.txt', package='chipenrich')
results = polyenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path,  method = 'polyenrich',
	locusdef = "nearest_tss", qc_plots = FALSE, out_name = NULL, n_cores = 1)
results.pe = results$results
print(results.pe[1:5,1:5])

## ---- fig.align='center', fig.cap='E2F4 peak distances to TSS', fig.height=6, fig.width=6, fig.show='hold', warning = FALSE, message = FALSE----
# Output in chipenrich and polyenrich
plot_dist_to_tss(peaks = peaks_E2F4, genome = 'hg19')

## ---- fig.align='center', fig.cap='E2F4 chipenrich spline without mappability', fig.height=6, fig.width=6, fig.show='hold', warning = FALSE, message = FALSE----
# Output in chipenrich
plot_chipenrich_spline(peaks = peaks_E2F4, locusdef = 'nearest_tss', genome = 'hg19')

## ---- fig.align='center', fig.cap='E2F4 polyenrich spline without mappability', fig.height=6, fig.width=6, fig.show='hold', warning = FALSE, message = FALSE----
# Output in chipenrich
plot_polyenrich_spline(peaks = peaks_E2F4, locusdef = 'nearest_tss', genome = 'hg19')

## ---- fig.align='center', fig.cap='H3K4me3 gene coverage', fig.height=6, fig.width=6, fig.show='hold', warning = FALSE, message = FALSE----
# Output in chipenrich
plot_gene_coverage(peaks = peaks_H3K4me3_GM12878, locusdef = 'nearest_tss',  genome = 'hg19')

## ------------------------------------------------------------------------
head(results$peaks)

## ------------------------------------------------------------------------
head(results$peaks_per_gene)

## ------------------------------------------------------------------------
head(results$results)

## ---- warning = FALSE, message = FALSE-----------------------------------
# Assessing if alpha = 0.05
gs_path = system.file('extdata','vignette_genesets.txt', package='chipenrich')
results = chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path,
	locusdef = "nearest_tss", qc_plots = FALSE, randomization = 'complete',
    out_name = NULL, n_cores = 1)
alpha = sum(results$results$P.value < 0.05) / nrow(results$results)
print(alpha)

