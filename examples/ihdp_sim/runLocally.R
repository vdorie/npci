numReps <- 250L

methods <- c("bart", "naive1", "naive2")
overlap <- c(TRUE, FALSE)
covariates <- c("full", "reduced")
verbose <- TRUE

baseResultInterval <- matrix(c(1L, numReps), 1L, 2L,
                             dimnames = list(NULL, c("start", "end")))

source("results.R")

for (i in seq_along(methods)) {
  for (j in seq_along(overlap)) {
    for (k in seq_along(covariates)) {
      existingResults <- getResultIntervals(methods[i], overlap[j], covariates[k])
      
      resultsIntervals <- intervalSubtraction(baseResultInterval, existingResults)
      if (nrow(resultsIntervals) == 0) next
      
      Sys.setenv(METHOD = methods[i])
      Sys.setenv(OVERLAP = if (overlap[j]) "true" else "false")
      Sys.setenv(COVARIATES = covariates[k])
      Sys.setenv(VERBOSE = verbose)
      
      for (l in seq_len(nrow(resultsIntervals))) {
        interval <- resultsIntervals[l,,drop = FALSE]
        
        Sys.setenv(START = interval[,"start"])
        Sys.setenv(END   = interval[,"end"])
        
        cat("running method '", methods[i], "' with overlap = ", overlap[j], " and covariates = '", covariates[k],
            "' from ", interval[,"start"], " to ", interval[,"end"], "\n", sep = "")
        
        source("runJob.R", local = new.env(parent = .GlobalEnv), keep.source = FALSE)
      }
    }
  }
}
