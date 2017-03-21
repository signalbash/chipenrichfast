context('Test read_bed()')

test_that('Uncompressed gff3 fly', {
	peak_path = system.file('extdata', 'test.gff3', package='chipenrich')
	peaks = read_bed(peak_path, genome = 'dm3')

	expect_equal(length(peaks), 19)
	expect_equal(unique(GenomeInfoDb::genome(peaks)), 'dm3')
})

test_that('Uncompressed gff3, no genome', {
	peak_path = system.file('extdata', 'test.gff3', package='chipenrich')
	peaks = read_bed(peak_path)

	expect_equal(length(peaks), 19)

})

test_that('Compressed gff3, no genome', {
	peak_path = system.file('extdata', 'test.gff3.gz', package='chipenrich')
	peaks = read_bed(peak_path)

	expect_equal(length(peaks), 19)
})

test_that('bedGraph human', {
	peak_path = system.file('extdata', 'test.bedGraph', package='chipenrich')
	peaks = read_bed(peak_path, genome = 'hg19')

	expect_equal(length(peaks), 4)
	expect_equal(unique(GenomeInfoDb::genome(peaks)), 'hg19')
})

test_that('BED3 with commented header', {
	peak_path = system.file('extdata', 'test_header.bed', package='chipenrich')
	peaks = suppressWarnings(read_bed(peak_path))

	expect_equal(length(peaks), 10)
	expect_warning(read_bed(peak_path), 'input regions overlap')
})

test_that('Uncompressed broadPeak', {
	peak_path = system.file('extdata', 'test.broadPeak', package='chipenrich')
	peaks = suppressWarnings(read_bed(peak_path))

	expect_equal(length(peaks), 27)
	expect_warning(read_bed(peak_path), 'input regions overlap')
})

test_that('Uncompressed narrowPeak', {
	peak_path = system.file('extdata', 'test.narrowPeak', package='chipenrich')
	peaks = read_bed(peak_path)

	expect_equal(length(peaks), 14)
})

test_that('Compressed broadPeak', {
	peak_path = system.file('extdata', 'test.broadPeak.gz', package='chipenrich')
	peaks = read_bed(peak_path, genome = 'hg19')

	expect_equal(length(peaks), 17)
})

test_that('Compressed narrowPeak', {
	peak_path = system.file('extdata', 'test.narrowPeak.gz', package='chipenrich')
	peaks = read_bed(peak_path)

	expect_equal(length(peaks), 11)
})

test_that('load_peaks() with and without genome', {
	peaks_df = data.frame(
		chr = c('chr1','chr2','chr3'),
		start = c(35,74,235),
		end = c(46,83,421),
		stringsAsFactors=F)

	# With genome
	peaks = load_peaks(peaks_df, genome = 'hg19')
	expect_equal(length(peaks), 3)

	# Without genome
	peaks = load_peaks(peaks_df)

	expect_equal(length(peaks), 3)
})

test_that('load_peaks() with extra columns', {
	peaks_df = data.frame(
		chr = c('chr1','chr2','chr3'),
		start = c(35,74,235),
		end = c(46,83,421),
		strand = c('+','-','+'),
		score = c(36, 747, 13),
		stringsAsFactors=F)
	peaks = load_peaks(peaks_df, keep.extra.columns = TRUE)

	expect_equal(length(peaks), 3)
})
