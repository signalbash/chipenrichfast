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

	expect_equal(min(sapply(as.list(filtered@set.gene), length)), 15)
	expect_equal(max(sapply(as.list(filtered@set.gene), length)), 100)
})

test_that('Filter genesets works for 5kb', {
	data(locusdef.hg19.5kb, package='chipenrich.data')
	data(geneset.GOBP.hsa, package='chipenrich.data')
	filtered = filter_genesets(geneset.GOBP.hsa, locusdef.hg19.5kb, 15, 100)

	expect_equal(min(sapply(as.list(filtered@set.gene), length)), 15)
	expect_equal(max(sapply(as.list(filtered@set.gene), length)), 100)
})

################################################################################
# Custom mappability definitions

test_that('Bad mappability header', {
	mappa_file = system.file('extdata', 'test_mappa_badhead.txt', package='chipenrich')
	expect_error(read_mappa(mappa_file), 'header must contain columns named')
})

test_that('Duplicate mappability gene_ids', {
	mappa_file = system.file('extdata', 'test_mappa_dup.txt', package='chipenrich')
	expect_error(read_mappa(mappa_file), 'duplicate gene_ids exist')
})

test_that('Negative mappability', {
	mappa_file = system.file('extdata', 'test_mappa_neg.txt', package='chipenrich')
	expect_error(read_mappa(mappa_file), 'mappability must be >= 0')
})

test_that('Too big mappability', {
	mappa_file = system.file('extdata', 'test_mappa_one.txt', package='chipenrich')
	expect_error(read_mappa(mappa_file), 'mappability must be <= 1')
})

test_that('Fine mappability', {
	mappa_file = system.file('extdata', 'test_mappa_good.txt', package='chipenrich')
	mappa = read_mappa(mappa_file)
	expect_true(all(mappa$gene_id == c(8487,84,91)))
})

################################################################################
# Custom locus definitions

ldef_noextras_file = system.file('extdata', 'test_ldef_noextras.txt', package='chipenrich')
ldef_geneid_file = system.file('extdata', 'test_ldef_geneid.txt', package='chipenrich')
ldef_symbol_file = system.file('extdata', 'test_ldef_symbol.txt', package='chipenrich')
ldef_nosymbol_file = system.file('extdata', 'test_ldef_nosymbol.txt', package='chipenrich')

test_that('Test no gene_id column throws error', {
	expect_error(setup_ldef(ldef_noextras_file), 'Custom locus definition must have column')
})

test_that('Test renaming geneid to gene_id and NA genome', {
	ldef = setup_ldef(ldef_geneid_file, genome = NA)

	expect_true(all(is.na(ldef@dframe$symbol)))
	expect_true(all(is.na(ldef@granges$symbol)))
})

test_that('Test renaming geneid to gene_id and hg19 genome', {
	ldef = setup_ldef(ldef_geneid_file, genome = 'hg19')

	expect_true(all(ldef@dframe$symbol == c('SAMD11','SAMD11','HES4','SAMD11','SAMD11')))
	expect_true(all(ldef@granges$symbol == c('SAMD11','SAMD11','HES4','SAMD11','SAMD11')))
})

test_that('Test nosymbol and NA genome', {
	ldef = setup_ldef(ldef_nosymbol_file, genome = NA)

	expect_true(all(is.na(ldef@dframe$symbol)))
	expect_true(all(is.na(ldef@granges$symbol)))
})

test_that('Test nosymbol and unsupported genome', {
	ldef = setup_ldef(ldef_nosymbol_file, genome = 'ab5')

	expect_true(all(is.na(ldef@dframe$symbol)))
	expect_true(all(is.na(ldef@granges$symbol)))
})

test_that('Test nosymbol and hg19 genome', {
	ldef = setup_ldef(ldef_nosymbol_file, genome = 'hg19')

	expect_true(all(ldef@dframe$symbol == c('A1BG','NAT2','CDH2')))
	expect_true(all(ldef@granges$symbol == c('A1BG','NAT2','CDH2')))
})

test_that('Test symbol and NA genome', {
	ldef = setup_ldef(ldef_symbol_file, genome = NA)

	expect_true(all(ldef@dframe$symbol == c('HI','BYE','SURE')))
	expect_true(all(ldef@granges$symbol == c('HI','BYE','SURE')))
})

test_that('Test symbol and hg19 genome', {
	ldef = setup_ldef(ldef_symbol_file, genome = 'hg19')

	expect_true(all(ldef@dframe$symbol == c('HI','BYE','SURE')))
	expect_true(all(ldef@granges$symbol == c('HI','BYE','SURE')))
})

test_that('Test symbol and unsupported genome', {
	ldef = setup_ldef(ldef_symbol_file, genome = 'ab5')

	expect_true(all(ldef@dframe$symbol == c('HI','BYE','SURE')))
	expect_true(all(ldef@granges$symbol == c('HI','BYE','SURE')))
})

################################################################################
# Custom genesets

# Test building the geneset itself
geneset_file = system.file('extdata', 'test_geneset.txt', package='chipenrich')
test_geneset = suppressMessages(setup_geneset(geneset_file))

test_that('Object of GeneSet class is returned',{
    expect_equal(class(test_geneset)[1], 'GeneSet')
})

test_that('Three genesets are created',{
    expect_equal(length(test_geneset@set.gene), 3)
})

test_that('Test names of genesets',{
    expect_equal(all(ls(test_geneset@set.gene) == c('GO:0035909','GO:0045822','GO:0045823')), TRUE)
})
