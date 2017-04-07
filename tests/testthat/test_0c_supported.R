context('Test supported functions')

test_that('Supported locus definitions works', {
	supported = supported_locusdefs()

	expect_equal(class(supported), 'data.frame')
})

test_that('Supported genesets works', {
	supported = supported_genesets()

	expect_equal(class(supported), 'data.frame')
})

test_that('Supported read lengths works', {
	supported = supported_read_lengths()

	expect_equal(class(supported), 'data.frame')
})

test_that('Supported genomes works', {
	supported = supported_genomes()

	expect_equal(class(supported), 'character')
	expect_equal(supported, c('danRer10', 'dm3', 'dm6', 'hg19', 'hg38', 'mm10', 'mm9', 'rn4', 'rn5', 'rn6'))
})

test_that('Supported methods works', {
	supported = supported_methods()

	expect_equal(class(supported), 'character')
	expect_equal(supported, c('chipenrich', 'chipenrich_fast', 'fet', 'broadenrich', 'polyenrich', 'polyenrich_fast', 'chipapprox'))
})
