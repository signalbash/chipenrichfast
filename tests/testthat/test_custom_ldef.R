context('Test User Supplied LocusDefinitions')

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
