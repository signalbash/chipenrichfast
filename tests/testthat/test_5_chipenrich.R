context('Test enrich functions')

test_that('Errors from no peaks', {
	ldef_file = system.file('extdata', 'test_ldef_symbol.txt', package = 'chipenrich')
	peak_file = system.file('extdata', 'test.broadPeak', package = 'chipenrich')
	gs_path = system.file('extdata','vignette_genesets.txt', package='chipenrich')

	expect_error(
		suppressWarnings(chipenrich(peaks = peak_file, genome = 'hg19', genesets = gs_path, locusdef = ldef_file, method = 'binomial', qc_plots = F, out_name = NULL, n_cores = 1)),
		'Intersection between peak midpoints and locus definition is empty')

	expect_error(
		suppressWarnings(chipenrich(peaks = peak_file, genome = 'hg19', genesets = gs_path, locusdef = ldef_file, qc_plots = F, out_name = NULL, n_cores = 1)),
		'Intersection between peak midpoints and locus definition is empty')

	expect_error(
		suppressWarnings(polyenrich(peaks = peak_file, genome = 'hg19', genesets = gs_path, locusdef = ldef_file, qc_plots = F, out_name = NULL, n_cores = 1)),
		'Intersection between peak midpoints and locus definition is empty')

	expect_error(
		suppressWarnings(broadenrich(peaks = peak_file, genome = 'hg19', genesets = gs_path, locusdef = ldef_file, qc_plots = F, out_name = NULL, n_cores=1)),
		'Intersection between peaks and locus definition is empty')
})

data(peaks_E2F4, package = 'chipenrich.data')
peaks_E2F4 = subset(peaks_E2F4, peaks_E2F4$chrom == 'chr1')

data(peaks_H3K4me3_GM12878, package = 'chipenrich.data')
peaks_H3K4me3_GM12878 = subset(peaks_H3K4me3_GM12878, peaks_H3K4me3_GM12878$chrom == 'chr1')

gs_path = system.file('extdata','vignette_genesets.txt', package='chipenrich')

test_that('Test chipenrich method without mappability', {
	results = suppressWarnings(chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path,
		locusdef = "nearest_tss", qc_plots = F, out_name = NULL, n_cores = 1))

	expect_equal(class(results), 'list')
})

test_that('Test chipenrich method with complete randomization', {
	results = suppressWarnings(chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path,
		locusdef = "nearest_tss", qc_plots = F, out_name = NULL, randomization = 'complete', n_cores = 1))

	expect_equal(class(results), 'list')
})

test_that('Test chipenrich method with bylength randomization', {
	results = suppressWarnings(chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path,
		locusdef = "nearest_tss", qc_plots = F, out_name = NULL, randomization = 'bylength', n_cores = 1))

	expect_equal(class(results), 'list')
})

test_that('Test chipenrich method with bylocation randomization', {
	results = suppressWarnings(chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path,
		locusdef = "nearest_tss", qc_plots = F, out_name = NULL, randomization = 'bylocation', n_cores = 1))

	expect_equal(class(results), 'list')
})

test_that('Test chipenrich method with mappability', {
	results = suppressWarnings(chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path,
		locusdef = "nearest_tss", mappability = '24', qc_plots = F,
		out_name = NULL,n_cores=1))

	expect_equal(class(results), 'list')
})

test_that('Test chipenrich_slow method', {
	results = suppressWarnings(chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path, method = 'chipenrich_slow',
		locusdef = "nearest_tss", qc_plots = F, out_name = NULL, n_cores = 1))

	expect_equal(class(results), 'list')
})

test_that('Test chipapprox method', {
    results = suppressWarnings(chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path, method = 'chipapprox',
    locusdef = "nearest_tss", qc_plots = F, out_name = NULL, n_cores = 1))
    
    expect_equal(class(results), 'list')
})

test_that('Test old chipapprox method', {
	results = suppressWarnings(chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path, method = 'chipapprox_old',
		locusdef = "nearest_tss", qc_plots = F, out_name = NULL, n_cores = 1))

	expect_equal(class(results), 'list')
})

test_that('Test polyenrich method', {
	results = suppressWarnings(polyenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path, method = 'polyenrich',
		locusdef = "nearest_tss", qc_plots = F, out_name = NULL, n_cores = 1))

	expect_equal(class(results), 'list')
})

test_that('Test polyenrich_slow method', {
	results = suppressWarnings(polyenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path, method = 'polyenrich_slow',
		locusdef = "nearest_tss", qc_plots = F, out_name = NULL, n_cores = 1))

	expect_equal(class(results), 'list')
})

test_that('Test polyapprox method', {
    results = suppressWarnings(polyenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path, method = 'polyapprox',
    locusdef = "nearest_tss", qc_plots = F, out_name = NULL, n_cores = 1))
    
    expect_equal(class(results), 'list')
})

test_that('Test broadenrich method', {
	results = suppressWarnings(broadenrich(peaks = peaks_H3K4me3_GM12878, genome = 'hg19', genesets = gs_path,
		locusdef = "nearest_tss", qc_plots = F, out_name = NULL, n_cores=1))

	expect_equal(class(results), 'list')
})

test_that('Test FET method', {
	results = suppressWarnings(chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path, locusdef = "5kb",
		method = "fet", fisher_alt = "two.sided", qc_plots = F, out_name = NULL))

	expect_equal(class(results), 'list')
})

test_that('Test binomial method', {
	results = suppressWarnings(chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path, locusdef = "5kb",
		method = "binomial", qc_plots = F, out_name = NULL))

	expect_equal(class(results), 'list')
})
