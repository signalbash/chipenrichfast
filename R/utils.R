# This is really a function to check for a valid method

# Recodes # of peaks to be:
# { 0, if no peaks
# { 1, if number of peaks is >= threshold argument
recode_peaks = function(num_peaks,threshold=1) {
  as.numeric(num_peaks >= threshold);
}

# Quick function to change column names of data frame
# (by name, not by position.)
change_names = function(dframe,name_list) {
  for (n in names(name_list)) {
    names(dframe)[names(dframe)==n] = name_list[[n]];
  }
  dframe;
}

# Remove a file extension from a file name.
strip_ext = function(file) {
  s = unlist(strsplit(file,".",fixed=T));
  s = s[1:(length(s)-1)];
  s = paste(s,collapse=".");
  return(s);
}

# Get a file's extension.
get_ext = function(file) {
  s = unlist(strsplit(file,".",fixed=T));
  return(tail(s,1));
}

# Checks to see if all elements of a list argument are possible.
check_arg = function(arg,possible_values,value=F) {
  if (!all(arg %in% possible_values)) {
    if (value) {
      return(arg[!arg %in% possible_values]);
    } else {
      return(F);
    }
  } else {
    return(T);
  }
}

genome_to_organism = function(genome) {
  code = substr(as.character(genome),1,2);
  if (code == 'mm') {
    org = 'mmu';
  } else if (code == 'hg') {
    org = 'hsa';
  } else if (code == 'rn') {
    org = 'rno';
  } else if (code == 'dm') {
  	org = 'dme';
  }
  else {
    org = NULL;
  }

  if (is.null(org)) {
    stop("Error: genome requested is not supported.");
  }

  org;
}

# Used to get eg2symbol mappings for custom locus definitions
genome_to_orgdb = function(genome = supported_genomes()) {
    genome = match.arg(genome)

    if(genome == 'hg19') {
        # Gives Entrez IDs
        egSYMBOL = org.Hs.eg.db::org.Hs.egSYMBOL
    } else if (genome == 'hg38') {
        # Gives Entrez IDs
        egSYMBOL = org.Hs.eg.db::org.Hs.egSYMBOL
    } else if (genome == 'mm9') {
        # Gives Entrez IDs
        egSYMBOL = org.Mm.eg.db::org.Mm.egSYMBOL
    } else if (genome == 'mm10') {
        # Gives Entrez IDs
        egSYMBOL = org.Mm.eg.db::org.Mm.egSYMBOL
    } else if (genome == 'rn4') {
        # Gives ENSEMBL IDs
        egSYMBOL = org.Rn.eg.db::org.Rn.egSYMBOL
    } else if (genome == 'rn5') {
        # Gives Entrez IDs
        egSYMBOL = org.Rn.eg.db::org.Rn.egSYMBOL
    } else if (genome == 'rn6') {
        # Gives Entrez IDs
        egSYMBOL = org.Rn.eg.db::org.Rn.egSYMBOL
    } else if (genome == 'dm3') {
        # Gives ENSEMBL IDs
        egSYMBOL = org.Dm.eg.db::org.Dm.egSYMBOL
    } else if (genome == 'dm6') {
        # Gives ENSEMBL IDs
        egSYMBOL = org.Dm.eg.db::org.Dm.egSYMBOL
    }

    ### Build Entrez ID to gene symbol mapping
        mapped_genes = AnnotationDbi::mappedkeys(egSYMBOL)
        eg2symbol = as.data.frame(egSYMBOL[mapped_genes])
        eg2symbol$gene_id = as.integer(eg2symbol$gene_id)

    return(eg2symbol)
}
