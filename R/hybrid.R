#' Running Hybrid test, either from scratch or using two results files
#'
#' Hybrid test is designed for people unsure of which test between ChIP-Enrich
#' and Poly-Enrich to use, so it takes information of both and gives adjusted
#' P-values. For more about ChIP- and Poly-Enrich, consult their corresponding
#' documentation.
#'
#' @section Hybrid p-values:
#' Given n tests that test for the same hypothesis, same Type I error rate, and
#' converted to p-values: \code{p_1, ..., p_n}, the Hybrid p-value is computed as:
#' \code{n*min(p_1, ..., p_n)}. This hybrid test will have at most the same
#' Type I error as any individual test, and if any of the tests have 100% power as
#' sample size goes to infinity, then so will the hybrid test.
#'
#' @section Function inputs:
#' Every input in hybridenrich is the same as in chipenrich and polyenrich. Inputs
#' unique to chipenrich are: num_peak_threshold; and inputs unique to polyenrich are:
#' weighting. Currently the test only supports running chipenrich and polyenrich, but
#' future plans will allow you to run any number of different support tests.
#'
#' @section Joining two results files
#' Combines two existing results files and returns one results file with hybrid
#' p-values and FDR included. Current allowed inputs are objects from any of
#' the supplied enrichment tests or a dataframe with at least the following columns:
#' \code{P.value, Geneset.ID}. Optional columns include: \code{Status}. Currently
#' we only allow for joining two results files, but future plans will allow you to join
#' any number of results files.


hybridenrich <- function(	peaks,
						out_name = "hybridenrich",
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
    if (!is.null(out_name)) {
        out_chip = sprintf("%s_chip",out_name)
        out_poly = sprintf("%s_poly",out_name)
    } else {
        out_chip = NULL
        out_poly = NULL
    }
    
	chip = chipenrich(
        peaks = peaks,
        out_name = out_chip,
        out_path = out_path,
        genome = genome,
        genesets = genesets,
        locusdef = locusdef,
        method = 'chipenrich',
        mappability = mappability,
        qc_plots = qc_plots,
        min_geneset_size = min_geneset_size,
        max_geneset_size = max_geneset_size,
        num_peak_threshold = num_peak_threshold,
        randomization = randomization,
        n_cores = n_cores)
        
    poly = polyenrich(
        peaks = peaks,
        out_name = out_poly,
        out_path = out_path,
        genome = genome,
        genesets = genesets,
        locusdef = locusdef,
        method = 'polyenrich',
        weighting = weighting,
        mappability = mappability,
        qc_plots = qc_plots,
        min_geneset_size = min_geneset_size,
        max_geneset_size = max_geneset_size,
        randomization = randomization,
        n_cores = n_cores)
        
    hybrid = hybrid.join(chip,poly)
    
    return(hybrid)
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
    if (!("Geneset.ID" %in% names(results1))) {
        stop("First object does not have Geneset.ID column")
    }
    if (!("Geneset.ID" %in% names(results2))) {
        stop("Second object does not have Geneset.ID column")
    }
    
    
    

    
    #Separate tree if the data does not have Status column
    if ("Status" %in% names(results1) & "Status" %in% names(results2)) {
        #Extract p-value and status of both tests
        Pvals1 = results1[,c("Geneset.ID","P.value","Status")]
        Pvals2 = results2[,c("Geneset.ID","P.value","Status")]
    
    } else {
        #Extract p-value only
        Pvals1 = results1[,c("Geneset.ID","P.value")]
        Pvals2 = results2[,c("Geneset.ID","P.value")]

    }

    #Merge by Geneset.ID
    PvalsH = merge(Pvals1, Pvals2, by="Geneset.ID")
    #If 0 remain, stop.
    if (nrow(PvalsH) == 0) {
        stop("No common genesets in the two datasets!")
    }
    message(sprintf("Total of %s common Geneset.IDs", nrow(PvalsH)))


    PvalsH$P.value.Hybrid = 2*pmin(PvalsH$P.value.x, PvalsH$P.value.y)
	
	#Run B-H to adjust for FDR for hybrid p-values
    PvalsH$FDR.Hybrid = stats::p.adjust(PvalsH$P.value.Hybrid, method = "BH")
    
    
    #Include enrich/depleted status if available and combine
    if ("Status" %in% names(results1) & "Status" %in% names(results2)) {
        PvalsH$Status.Hybrid = ifelse(PvalsH$Status.x == PvalsH$Status.y, PvalsH$Status.x, "Inconsistent")
        #Combine both results files together and append hybrid p-value and FDR
        resultsH = merge(results1[,-which(colnames(results1) %in% c("P.value","Status"))], PvalsH, by = "Geneset.ID")
    } else {
        resultsH = merge(results1[,-which(colnames(results1) %in% c("P.value"))], PvalsH, by = "Geneset.ID")
    }
    
    
    #Reorder columns?????
	
	#Output final results
	return(resultsH)
	
}