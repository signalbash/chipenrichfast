context('Test setup functions')

################################################################################
# Filter genesets

test_that('Filter genesets errors', {
	data(locusdef.hg19.nearest_tss, package='chipenrich.data')
	expect_error(filter_genesets('test', locusdef.hg19.nearest_tss, 20, 100), 'gs_obj not of class GeneSet')
})

test_that('Filter locus def errors', {
	data(geneset.GOBP.hsa, package='chipenrich.data')
	expect_error(filter_genesets(geneset.GOBP.hsa, 'test', 20, 100), 'ldef_obj not of class LocusDefinition')
})

test_that('Filter genesets works for nearest_tss', {
	data(locusdef.hg19.nearest_tss, package='chipenrich.data')
	data(geneset.GOBP.hsa, package='chipenrich.data')
	filtered = filter_genesets(geneset.GOBP.hsa, locusdef.hg19.nearest_tss, 15, 100)

	expect_equal(min(sapply(as.list(filtered@set.gene), length)),15)
	expect_equal(max(sapply(as.list(filtered@set.gene), length)), 100)
})

test_that('Filter genesets works for 5kb', {
	data(locusdef.hg19.5kb, package='chipenrich.data')
	data(geneset.GOBP.hsa, package='chipenrich.data')
	filtered = filter_genesets(geneset.GOBP.hsa, locusdef.hg19.5kb, 20, 100)

	expect_equal(min(sapply(as.list(filtered@set.gene), length)), 20)
	expect_equal(max(sapply(as.list(filtered@set.gene), length)), 100)
})
