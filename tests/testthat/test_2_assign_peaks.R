context('Test peak assignment functions')

################################################################################
# Test human locus definitions

# Load the required data for the test
hg19_ldefs = grep('locusdef.hg19', data(package = 'chipenrich.data')$results[,3], value=T)
data(list = hg19_ldefs, package = 'chipenrich.data')
tss = get(data('tss.hg19', package = 'chipenrich.data'))

file = system.file('extdata', 'test_assign.bed', package = 'chipenrich')
peaks = suppressMessages(read_bed(file))

# Expected values for peak assignments
expected_genes = list(
	'locusdef.hg19.10kb' = c(100287102),
	'locusdef.hg19.10kb_outside' = c(641702),
	'locusdef.hg19.10kb_outside_upstream' = c(641702),
	'locusdef.hg19.1kb' = c(100287102),
	'locusdef.hg19.1kb_outside' = c(100287102, 641702),
	'locusdef.hg19.1kb_outside_upstream' = c(100287102, 641702),
	'locusdef.hg19.5kb' = c(100287102),
	'locusdef.hg19.5kb_outside' = c(641702),
	'locusdef.hg19.5kb_outside_upstream' = c(641702),
	'locusdef.hg19.exon' = c(100287102),
	'locusdef.hg19.intron' = c(100287102),
	'locusdef.hg19.nearest_gene' = c(100287102, 641702),
	'locusdef.hg19.nearest_tss' = c(100287102, 641702)
)

# Do the tests in a loop for each locus definition of hg19
for(ldef in hg19_ldefs) {
	test_assign = assign_peaks(peaks = peaks, locusdef = get(ldef), tss = tss)

	expect_true(base::setequal(test_assign$gene_id, expected_genes[[ldef]]),
		info = sprintf('Test assign_peaks(): %s', ldef))
}

expected_overlaps = list(
	'locusdef.hg19.10kb' = c(500,100,100,80),
	'locusdef.hg19.10kb_outside' = c(921),
	'locusdef.hg19.10kb_outside_upstream' = c(921),
	'locusdef.hg19.1kb' = c(127,100,100),
	'locusdef.hg19.1kb_outside' = c(374,1000),
	'locusdef.hg19.1kb_outside_upstream' = c(374,1000),
	'locusdef.hg19.5kb' = c(500,100,100),
	'locusdef.hg19.5kb_outside' = c(1000,1000),
	'locusdef.hg19.5kb_outside_upstream' = c(1000,1000),
	'locusdef.hg19.exon' = c(100),
	'locusdef.hg19.intron' = c(100),
	'locusdef.hg19.nearest_gene' = c(500,100,100,1000),
	'locusdef.hg19.nearest_tss' = c(500,100,100,1000)
)

# Do the tests in a loop for each locus definition of hg19
for(ldef in hg19_ldefs) {
	test_assign = assign_peak_segments(peaks = peaks, locusdef = get(ldef))

	expect_true(all(test_assign$peak_overlap == expected_overlaps[[ldef]]),
		info = sprintf('Test assign_peak_segments(): %s', ldef))
}
