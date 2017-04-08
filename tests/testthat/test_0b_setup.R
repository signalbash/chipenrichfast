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

test_that('Errors in setup_locusdef()', {
	expect_error(setup_locusdef('blue', genome = 'hg19'), 'invalid genome / definition combination')
	expect_error(setup_locusdef('nearest_tss', genome = 'hg19', randomization = 'blue'), 'Invalid randomization')
})

test_that('Errors in setup_genesets()', {
	expect_error(setup_genesets(gs_codes = 'GOBP', ldef_obj = 'blue', genome = 'hg19', min_geneset_size = 15, max_geneset_size = 2000), 'ldef_obj not of class LocusDefinition')

	ldef = setup_locusdef('nearest_tss', genome = 'hg19')$ldef
	expect_error(setup_genesets(gs_codes = 'blue', ldef_obj = ldef, genome = 'hg19', min_geneset_size = 15, max_geneset_size = 2000), 'Invalid organism / geneset combination requested')
})

test_that('Errors in setup_mappa()', {
	mappa_file = system.file('extdata', 'test_mappa_good.txt', package = 'chipenrich')
	ldef_file = system.file('extdata', 'test_ldef_symbol.txt', package = 'chipenrich')
	
	ldef = setup_locusdef('nearest_tss', genome = 'hg19')$ldef
	custom_ldef = setup_locusdef(ldef_file, genome = 'hg19')$ldef

	expect_error(setup_mappa(mappa_code = '24', genome = 'hg19', ldef_code = 'blue', ldef_obj = 'blue'), 'ldef_obj not of class LocusDefinition')
	expect_error(setup_mappa(mappa_code = mappa_file, genome = 'hg19', ldef_code = 'nearest_tss', ldef_obj = ldef), 'your mappability genes and locus definition genes overlap')
	expect_error(setup_mappa(mappa_code = '24', genome = 'hg19', ldef_code = ldef_file, ldef_obj = custom_ldef), 'Built-in mappability cannot be used with a user-defined locus definition')
})
