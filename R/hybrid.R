#User gives method names, and then all of them are run
hybrid.enrich <- function(	peaks,
						out_name = "polyenrich",
						out_path = getwd(),
						genome = supported_genomes(),
						genesets = c(
							'GOBP',
							'GOCC',
							'GOMF'),
						locusdef = "nearest_tss",
						methods = c('chipenrich','polyenrich'),
						weighting = NULL,
						mappability = NULL,
						qc_plots = TRUE,
						min_geneset_size = 15,
						max_geneset_size = 2000,
						num_peak_threshold = 1,
						randomization = NULL,
						n_cores = 1
) {
	NULL
}



#User gives results objects (object or just results), and appends hybrid results
hybrid.join <- function(test1, test2) {
	#Check if they inputed the test object or just the results file, checked by seeing
    # if object has a $results part
    if ("results" %in% names(test1)) {
        #If entire object, extract the results section
        results1 = test1$results
    } else if ("P.value" %in% names(test1)) {
        results1 = test1
    } else {
        stop("First object is not a valid output or does not have P.value column")
    }
	
    if ("results" %in% names(test2)) {
        #If entire object, extract the results section
        results2 = test2$results
    } else if ("P.value" %in% names(test2)) {
        results2 = test2
    } else {
        stop("Second object is not a valid output or does not have P.value column")
    }
    
    #Check for Geneset.ID column
    
    #Separate tree if the data does not have Status column
        
	#Extract p-value of both tests and calculate hybrid
	Pvals1 = results1[,c("Geneset.ID","P.value","Status")]
    Pvals2 = results2[,c("Geneset.ID","P.value","Status")]
    
    PvalsH = merge(Pvals1, Pvals2, by="Geneset.ID")
    
    PvalsH$P.value.Hybrid = 2*pmin(PvalsH$P.value.x, PvalsH$P.value.y)
	
	#Run B-H to adjust for FDR for hybrid p-values
    PvalsH$FDR.Hybrid = stats::p.adjust(PvalsH$P.value.Hybrid, method = "BH")
    
    #Include enrich/depleted status
    PvalsH$Status.Hybrid = ifelse(PvalsH$Status.x == PvalsH$Status.y, PvalsH$Status.x, "Inconsistent")
	
	#Combine both results files together and append hybrid p-value and FDR
	merge(results1)
	
	#Output final results
	return(out)
	
}