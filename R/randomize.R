# Randomize the locus definition @dframe and rebuild the @granges and @chrom2iranges
randomize_ldef_bylocation = function(ldef, resolution=50) {
	message('Randomizing locus definition in location bins...')
	# Extract GRanges representation
	ldef_gr = ldef@granges

	# Split by chromosome
	ldef_grl = S4Vectors::split(ldef_gr, GenomeInfoDb::seqnames(ldef_gr))

	# Within each chromosome:
	# 1. Form groups based on the number of genes on that chromosome
	# 2. Split the chromosome GRanges by the groups (GRangesList)
	# 3. Within each group within each chromosome, scramble the mcols
	# 4. Collapse the groups so the result is a GRanges for the chromosome again (empty GRanges for chromosomes without genes).
	re_ldef_grl = S4Vectors::endoapply(ldef_grl, function(gr) {
		if(length(gr) > 0) {
			# Create groups within the chromosome
			group = floor( (seq_along(gr) + (resolution - 1)) / resolution )
			# Split by the group
			grl = S4Vectors::splitAsList(gr, group)
			# Rearrange mcols within each group
			re_grl = S4Vectors::endoapply(grl, function(grg){
				GenomicRanges::mcols(grg) = GenomicRanges::mcols(grg)[sample(seq_along(grg), length(grg)),]
				return(grg)
			})
			# Collapse the GRangesList into the scrambled GRanges for the chromosome
			re_gr = BiocGenerics::unlist(re_grl, use.names = FALSE)
		} else {
			# Need this for chromosomes without genes
			re_gr = GenomicRanges::GRanges()
		}

		return(re_gr)
	})
	# Collapse the GRangesList across the chromosomes
	re_ldef_gr = unlist(re_ldef_grl, use.names = FALSE)
	# Make sure it's sorted
	re_ldef_gr = sort(re_ldef_gr)

	# Construct data.frame for new locus definition
	re_ldef_df = data.frame(re_ldef_gr, stringsAsFactors = FALSE)
	re_ldef_df = re_ldef_df[, c('seqnames','start','end','gene_id','symbol')]
	colnames(re_ldef_df) = c('chr','start','end','gene_id','symbol')

	# Reassign the dframe and granges to the ldef
	ldef@dframe = re_ldef_df
	ldef@granges = re_ldef_gr

	return(ldef)
}

# Randomize ppg after all additions have been made across all genes
randomize_ppg_all = function(ppg) {
	ppg = ppg[order(ppg$length),]
	rownames(ppg) = 1:nrow(ppg)

	reordering = sample(1:nrow(ppg), nrow(ppg))
	ppg = data.frame('gene_id'=ppg$gene_id, ppg[reordering,2:ncol(ppg)], stringsAsFactors = FALSE)

	return(ppg)
}

# Randomize ppg after all additions have been made across all genes
randomize_ldef_complete = function(ldef) {
	message('Randomizing locus definition...')
	gr = ldef@granges
	df = ldef@dframe

	shuffle = base::sample(seq_along(gr), length(gr))

	GenomicRanges::mcols(gr) = GenomicRanges::mcols(gr)[shuffle, ]
	df[, c('gene_id', 'symbol')] = df[shuffle, c('gene_id', 'symbol')]

	ldef@granges = gr
	ldef@dframe = df

	return(ldef)
}

# Randomize ppg after all additions have been made within length bins
randomize_ppg_length = function(ppg) {
	ppg = ppg[sample(1:nrow(ppg),nrow(ppg)),]
	ppg = ppg[order(ppg$length),]
	rownames(ppg) = 1:nrow(ppg)

	group = floor(as.numeric(rownames(ppg))+99)/100
	group = floor(group)

	split_ppg = split(ppg, group)
	split_ppg = lapply(split_ppg, function(bin){
		reordering = sample(1:nrow(bin), nrow(bin))

		data.frame('gene_id'=bin$gene_id, bin[reordering,2:ncol(bin)], stringsAsFactors = FALSE)
	})
	ppg = Reduce(rbind, split_ppg)

	return(ppg)
}

randomize_ldef_bylength = function(ldef, resolution = 100) {
	message('Randomizing locus definition in locus length bins...')
	gr = ldef@granges
	df = ldef@dframe

	widths = data.frame(
		original_idx = seq_along(gr),
		width = GenomicRanges::width(gr),
		stringsAsFactors = FALSE)
	widths = widths[order(widths$width), ]
	rownames(widths) = seq_along(gr)

	group = floor(as.numeric(rownames(widths))+ (resolution - 1)) / resolution
	group = floor(group)

	split_widths = split(widths, group)
	split_widths = lapply(split_widths, function(bin){
		shuffle = base::sample(1:nrow(bin), nrow(bin))
		bin$shuffle_idx = bin$original_idx[shuffle]

		return(bin)
	})
	widths = Reduce(rbind, split_widths)

	GenomicRanges::mcols(gr)[widths$original_idx, ] = GenomicRanges::mcols(gr)[widths$shuffle_idx, ]
	df[widths$original_idx, c('gene_id', 'symbol')] = df[widths$shuffle_idx, c('gene_id', 'symbol')]

	ldef@granges = gr
	ldef@dframe = df

	return(ldef)
}
