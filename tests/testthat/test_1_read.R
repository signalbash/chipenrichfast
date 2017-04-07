context('Test read functions')

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
	expect_error(read_ldef(ldef_noextras_file), 'Custom locus definition must have column')
})

test_that('Test renaming geneid to gene_id and NA genome', {
	ldef = read_ldef(ldef_geneid_file, genome = NA)

	expect_true(all(is.na(ldef@dframe$symbol)))
	expect_true(all(is.na(ldef@granges$symbol)))
})

test_that('Test renaming geneid to gene_id and hg19 genome', {
	ldef = read_ldef(ldef_geneid_file, genome = 'hg19')

	expect_true(all(ldef@dframe$symbol == c('SAMD11','SAMD11','HES4','SAMD11','SAMD11')))
	expect_true(all(ldef@granges$symbol == c('SAMD11','SAMD11','HES4','SAMD11','SAMD11')))
})

test_that('Test nosymbol and NA genome', {
	ldef = read_ldef(ldef_nosymbol_file, genome = NA)

	expect_true(all(is.na(ldef@dframe$symbol)))
	expect_true(all(is.na(ldef@granges$symbol)))
})

test_that('Test nosymbol and unsupported genome', {
	ldef = read_ldef(ldef_nosymbol_file, genome = 'ab5')

	expect_true(all(is.na(ldef@dframe$symbol)))
	expect_true(all(is.na(ldef@granges$symbol)))
})

test_that('Test nosymbol and hg19 genome', {
	ldef = read_ldef(ldef_nosymbol_file, genome = 'hg19')

	expect_true(all(ldef@dframe$symbol == c('A1BG','NAT2','CDH2')))
	expect_true(all(ldef@granges$symbol == c('A1BG','NAT2','CDH2')))
})

test_that('Test symbol and NA genome', {
	ldef = read_ldef(ldef_symbol_file, genome = NA)

	expect_true(all(ldef@dframe$symbol == c('HI','BYE','SURE')))
	expect_true(all(ldef@granges$symbol == c('HI','BYE','SURE')))
})

test_that('Test symbol and hg19 genome', {
	ldef = read_ldef(ldef_symbol_file, genome = 'hg19')

	expect_true(all(ldef@dframe$symbol == c('HI','BYE','SURE')))
	expect_true(all(ldef@granges$symbol == c('HI','BYE','SURE')))
})

test_that('Test symbol and unsupported genome', {
	ldef = read_ldef(ldef_symbol_file, genome = 'ab5')

	expect_true(all(ldef@dframe$symbol == c('HI','BYE','SURE')))
	expect_true(all(ldef@granges$symbol == c('HI','BYE','SURE')))
})

################################################################################
# Custom genesets

# Test building the geneset itself
geneset_file = system.file('extdata', 'test_geneset.txt', package='chipenrich')
test_geneset = suppressMessages(read_geneset(geneset_file))

test_that('Object of GeneSet class is returned',{
    expect_equal(class(test_geneset)[1], 'GeneSet')
})

test_that('Three genesets are created',{
    expect_equal(length(test_geneset@set.gene), 3)
})

test_that('Test names of genesets',{
    expect_equal(all(ls(test_geneset@set.gene) == c('GO:0035909','GO:0045822','GO:0045823')), TRUE)
})

################################################################################
# read_bed()

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

test_that('Text input header', {
	peak_path = system.file('extdata', 'test_header.txt', package='chipenrich')
	peaks = suppressWarnings(read_bed(peak_path))

	expect_equal(length(peaks), 10)
	expect_warning(read_bed(peak_path), 'input regions overlap')
})

test_that('Text input no header', {
	peak_path = system.file('extdata', 'test_noheader.txt', package='chipenrich')
	peaks = suppressWarnings(read_bed(peak_path))

	expect_equal(length(peaks), 10)
	expect_warning(read_bed(peak_path), 'input regions overlap')
})

################################################################################
# load_peaks()

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
