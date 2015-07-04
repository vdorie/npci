numReps <- 250L

methods <- c("bart", "naive1", "naive2")

baseResultInterval <- matrix(c(1L, totalNumReps), 1L, 2L,
                             dimnames = list(NULL, c("start", "end")))

source("results.R")

for (method in methods) {
  existingResults <- getResultIntervals(method)
  
  resultsIntervals <- intervalSubtraction(baseResultInterval, existingResults)
  if (nrow(resultsIntervals) == 0) next
  
  Sys.setenv(METHOD = method)
  
  for (j in seq_len(nrow(resultsIntervals))) {
    interval <- resultsIntervals[j,, drop = FALSE]
    
    Sys.setenv(START = interval[, "start"])
    Sys.setenv(END   = interval[, "end"])
    
    cat("running method '", method, "' from ", interval[, "start"], " to ", interval[, "end"], "\n", sep = "")
    
    source("runJob.R", local = TRUE, keep.source = FALSE)
  }
}
