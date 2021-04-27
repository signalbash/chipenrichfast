context('Test peak assignment functions')

################################################################################
# Test no assigned peaks

test_that('No peaks assigned', {

	ldef_file = system.file('extdata', 'test_ldef_symbol.txt', package = 'chipenrich')
	peak_file = system.file('extdata', 'test.broadPeak', package = 'chipenrich')

	peaks = suppressWarnings(read_bed(peak_file))
	ldef_list = setup_locusdef(ldef_file, genome = 'hg19')

	expect_error(assign_peaks(peaks, ldef_list$ldef, ldef_list$tss), 'Intersection between peak midpoints and locus definition is empty')
	expect_error(assign_peak_segments(peaks, ldef_list$ldef), 'Intersection between peaks and locus definition is empty')
})

################################################################################
# Test human locus definitions
