#User gives method names, and then all of them are run
hybrid.run <- function(	peaks,
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
	#Check if they inputed the test object or just the results file
	NULL
	#If entire object, extract the results section
	#    object$results
	
	#Extract p-value of both tests
	NULL
	
	#Run B-H to adjust for FDR for hybrid p-values
	NULL
	
	#Combine both results files together and append hybrid p-value and FDR
	NULL
	
	#Output final results
	NULL
	
}