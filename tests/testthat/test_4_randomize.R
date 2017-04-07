context('Test randomize functions')

################################################################################
# Test complete randomization
test_that('Test complete locusdef randomization', {
	data(locusdef.dm6.nearest_tss, package = 'chipenrich.data')

	randomized = randomize_ldef_complete(locusdef.dm6.nearest_tss)

	expect_equal(nrow(locusdef.dm6.nearest_tss@dframe), nrow(randomized@dframe))
	expect_equal(length(locusdef.dm6.nearest_tss@granges), length(randomized@granges))
	expect_equal(locusdef.dm6.nearest_tss@genome.build, randomized@genome.build)
	expect_equal(locusdef.dm6.nearest_tss@organism, randomized@organism)
	expect_true(all(locusdef.dm6.nearest_tss@granges$gene_id %in% randomized@granges$gene_id))
})

################################################################################
# Test bylength randomization by location
test_that('Test locusdef randomization by length', {
	data(locusdef.dm6.nearest_tss, package = 'chipenrich.data')

	randomized = randomize_ldef_bylength(locusdef.dm6.nearest_tss, resolution = 100)

	expect_equal(nrow(locusdef.dm6.nearest_tss@dframe), nrow(randomized@dframe))
	expect_equal(length(locusdef.dm6.nearest_tss@granges), length(randomized@granges))
	expect_equal(locusdef.dm6.nearest_tss@genome.build, randomized@genome.build)
	expect_equal(locusdef.dm6.nearest_tss@organism, randomized@organism)
	expect_true(all(locusdef.dm6.nearest_tss@granges$gene_id %in% randomized@granges$gene_id))
})

################################################################################
# Test locusdef randomization by location
test_that('Test locusdef randomization by location', {
	data(locusdef.dm6.nearest_tss, package = 'chipenrich.data')

	randomized = randomize_ldef_bylocation(locusdef.dm6.nearest_tss, resolution = 100)

	expect_equal(nrow(locusdef.dm6.nearest_tss@dframe), nrow(randomized@dframe))
	expect_equal(length(locusdef.dm6.nearest_tss@granges), length(randomized@granges))
	expect_equal(locusdef.dm6.nearest_tss@genome.build, randomized@genome.build)
	expect_equal(locusdef.dm6.nearest_tss@organism, randomized@organism)
	expect_true(all(locusdef.dm6.nearest_tss@granges$gene_id %in% randomized@granges$gene_id))
})

################################################################################
# Test randomizations without binning
test_that('Test ppg randomization without length bins', {
	ppg_file = system.file('extdata', 'test_ppg.txt', package = 'chipenrich')
	ppg = read.table(ppg_file, header = T, sep='\t', stringsAsFactors=F)

	randomized = randomize_ppg_all(ppg)

	expect_equal(mean(ppg$length), mean(randomized$length))
	expect_equal(mean(ppg$log10_length), mean(randomized$log10_length))
	expect_equal(mean(ppg$num_peaks), mean(randomized$num_peaks))
	expect_true(all(ppg$gene_id %in% randomized$gene_id))
})

################################################################################
# Test randomization in bins
test_that('Test ppg randomization within length bins', {
	ppg_file = system.file('extdata', 'test_ppg.txt', package = 'chipenrich')
	ppg = read.table(ppg_file, header = T, sep='\t', stringsAsFactors=F)

	randomized = randomize_ppg_length(ppg)

	expect_equal(mean(ppg$length), mean(randomized$length))
	expect_equal(mean(ppg$log10_length), mean(randomized$log10_length))
	expect_equal(mean(ppg$num_peaks), mean(randomized$num_peaks))
	expect_true(all(ppg$gene_id %in% randomized$gene_id))
})
