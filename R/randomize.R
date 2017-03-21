# Randomize the locus definition @dframe and rebuild the @granges and @chrom2iranges
randomize_locusdef = function(ldef, resolution=50) {
	### TODO: Remove the reducing code once we transition to new chipenrich.data

	message('Extracting definition data frame..')
	ldef_df = ldef@dframe

	message('Splitting data frame on chrom..')
	ldef_chrom = split(ldef_df, f=ldef_df$chrom)

	message('Working within each chrom..')
	ldef_chrom = lapply(ldef_chrom, function(lc) {
		chrom = unique(lc$chrom)

		message(paste('	On chrom', chrom))

		# Collapse consecutive rows of the ldef when the gene_id is the same
			# Split by gene_id
			message('		Splitting chrom on gene id..')
			lc_gene_id = split(lc, lc$gene_id)

			# For each gene_id df, create an IRanges object to apply IRanges::reduce to
			message('		Applying IRanges::reduce to each gene id group..')
			lc_gene_id_ir = lapply(lc_gene_id, function(lcg) {
				gene_id = unique(lcg$gene_id)

				ir = IRanges::IRanges(start=lcg$start, end = lcg$end, names=lcg$gene_id)
				ir = IRanges::reduce(ir)
				names(ir) = rep.int(gene_id, length(ir))

				return(as.data.frame(ir, stringsAsFactors=F))
			})

		# Now put all the data.framed IRange objects back together
		message('		Collapsing reduced definition..')
		lc_collapsed = Reduce(rbind, lc_gene_id_ir)

		# Add chromosome column back, rename names column, and sort columns
		message('		Formatting reduced definition..')
		lc_collapsed$chrom = chrom
		colnames(lc_collapsed) = c('start','end','width','gene_id','chrom')
		lc_collapsed = lc_collapsed[,c('gene_id','chrom','start','end')]

		# Sort lc_collapsed by the starting position and rename rownames
		message('		Sorting reduced definition..')
		lc_collapsed = lc_collapsed[order(lc_collapsed$start),]
		rownames(lc_collapsed) = 1:nrow(lc_collapsed)

		# Form groups
		group = floor(as.numeric(rownames(lc_collapsed))+(resolution-1))/resolution
		group = floor(group)

		# Split the chromosome into parts by group
		message(paste('		Shuffling within bins of', resolution, 'genes'))
		split_lc = split(lc_collapsed, group)
		split_lc = lapply(split_lc, function(bin){
			reordering = sample(1:nrow(bin), nrow(bin))

			# Scramble gene_ids
			data.frame('gene_id' = bin$gene_id[reordering], bin[,2:ncol(bin)], stringsAsFactors=F)
		})
		lc = Reduce(rbind, split_lc)

		return(lc)
	})

	message('Done rejiggering genomic locations for data frame..')
	ldef_df = Reduce(rbind, ldef_chrom)

	message('Creating new GenomicRanges object..')
	ldef_gr = GenomicRanges::GRanges(
		seqnames = ldef_df$chrom,
		ranges = IRanges::IRanges(start=ldef_df$start, end=ldef_df$end),
		names = ldef_df$gene_id
	)

	message('Creating new IRanges object..')
	chroms = c(paste('chr',1:22,sep=''),'chrX','chrY')

	chr_list = chroms
	names(chr_list) = chroms

	ldef_ir = lapply(chr_list, function(chr) {
		sub_ldef_df = subset(ldef_df, ldef_df$chrom==chr)

		ir = IRanges::IRanges(start=sub_ldef_df$start, end=sub_ldef_df$end, names=sub_ldef_df$gene_id)
		return(ir)
	})

	ldef@dframe = ldef_df
	ldef@granges = ldef_gr
	ldef@chrom2iranges = ldef_ir

	return(ldef)
}

# Randomize ppg after all additions have been made across all genes
randomize_ppg_all = function(ppg) {
	ppg = ppg[order(ppg$length),]
	rownames(ppg) = 1:nrow(ppg)

	reordering = sample(1:nrow(ppg), nrow(ppg))
	ppg = data.frame('gene_id'=ppg$gene_id, ppg[reordering,2:ncol(ppg)], stringsAsFactors=F)

	return(ppg)
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

		data.frame('gene_id'=bin$gene_id, bin[reordering,2:ncol(bin)], stringsAsFactors=F)
	})
	ppg = Reduce(rbind, split_ppg)

	return(ppg)
}
