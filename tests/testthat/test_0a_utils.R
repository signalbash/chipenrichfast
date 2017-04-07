context('Test utility functions')

test_that('Recode peaks', {
    num_peaks = c(1,4,7,1,8,0,5)

    expect_true(all(recode_peaks(num_peaks) == c(1,1,1,1,1,0,1)))
    expect_true(all(recode_peaks(num_peaks, threshold = 2) == c(0,1,1,0,1,0,1)))
})

test_that('Genome to organism', {
    expect_error(genome_to_organism('test'), "'arg' should be one of")
    expect_equal(genome_to_organism('dm3'), 'dme')
    expect_equal(genome_to_organism('dm6'), 'dme')
    expect_equal(genome_to_organism('danRer10'), 'dre')
    expect_equal(genome_to_organism('hg19'), 'hsa')
    expect_equal(genome_to_organism('hg38'), 'hsa')
    expect_equal(genome_to_organism('mm9'), 'mmu')
    expect_equal(genome_to_organism('mm10'), 'mmu')
    expect_equal(genome_to_organism('rn4'), 'rno')
    expect_equal(genome_to_organism('rn5'), 'rno')
    expect_equal(genome_to_organism('rn6'), 'rno')
})

test_that('Genome to eg2symbol', {
    expect_error(genome_to_orgdb('test'), "'arg' should be one of")
    expect_equal(class(genome_to_orgdb('dm6')), 'data.frame')
    expect_equal(class(genome_to_orgdb('danRer10')), 'data.frame')
    expect_equal(class(genome_to_orgdb('hg19')), 'data.frame')
    expect_equal(class(genome_to_orgdb('mm10')), 'data.frame')
    expect_equal(class(genome_to_orgdb('rn6')), 'data.frame')
})
