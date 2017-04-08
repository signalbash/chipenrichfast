context('Test QC plot functions')

data(peaks_E2F4, package = 'chipenrich.data')
data(peaks_H3K4me3_GM12878, package = 'chipenrich.data')

test_that('Distance to TSS plot works', {
	plot = suppressWarnings(plot_dist_to_tss(peaks = peaks_E2F4, genome = 'hg19'))

	expect_equal(class(plot), 'trellis')
})

test_that('ChIP-Enrich spline plot works', {
	plot = suppressWarnings(plot_chipenrich_spline(peaks = peaks_E2F4, locusdef = 'nearest_tss', genome = 'hg19'))

	expect_equal(class(plot), 'trellis')
})

test_that('ChIP-Enrich spline plot with mappability works', {
	plot = suppressWarnings(plot_chipenrich_spline(peaks = peaks_E2F4, locusdef = 'nearest_tss', genome = 'hg19', mappability = '50'))

	expect_equal(class(plot), 'trellis')
})

test_that('Poly-Enrich spline plot works', {
	plot = suppressWarnings(plot_polyenrich_spline(peaks = peaks_E2F4, locusdef = 'nearest_tss', genome = 'hg19'))

	expect_equal(class(plot), 'trellis')
})

test_that('Broad-Enrich coverage plot works', {
	plot = suppressWarnings(plot_gene_coverage(peaks = peaks_H3K4me3_GM12878, locusdef = 'nearest_tss',  genome = 'hg19'))

	expect_equal(class(plot), 'trellis')
})
