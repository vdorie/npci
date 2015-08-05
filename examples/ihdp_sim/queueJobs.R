## for use with TORQUE/PBS/qsub

if (require(npci, quietly = TRUE) == FALSE)
  stop("npci package not available")

totalNumReps <- 250L

#methods <- c("bart", "naive1", "naive2")
#numRepsPerProcess <- c(250L, 5L, 10L)
#overlap <- c(TRUE, FALSE)
#covariates <- c("full", "reduced")

methods <- "grouped"
numRepsPerProcess <- 5L
overlap <- TRUE
covariates <- "full"

verbose <- TRUE

if (!dir.exists("jobs")) dir.create("jobs")

baseResultInterval <- matrix(c(1L, totalNumReps), 1L, 2L,
                             dimnames = list(NULL, c("start", "end")))

source("results.R")

jobIter <- 1L

for (i in seq_along(methods)) {
  for (j in seq_along(overlap)) {
    for (k in seq_along(covariates)) {
      prefix <- getPrefix(methods[i], overlap[j], covariates[k])
  
      existingResults <- getResultIntervals(methods[i], overlap[j], covariates[k])
      
      resultsIntervals <- intervalSubtraction(baseResultInterval, existingResults)
      if (nrow(resultsIntervals) == 0) next
      
      for (l in seq_len(nrow(resultsIntervals))) {
        interval <- resultsIntervals[l,,drop = FALSE]
        numRepsPerInterval <- unname(interval[,"end"] - interval[,"start"] + 1L)
        
        numProcesses <- numRepsPerInterval %/% numRepsPerProcess[i] +
          if (numRepsPerInterval %% numRepsPerProcess[i] != 0) 1L else 0L
        numReps <- rep(numRepsPerInterval %/% numProcesses, numProcesses) +
          c(rep(1L, numRepsPerInterval %% numProcesses), rep(0L, numProcesses - numRepsPerInterval %% numProcesses))
        
        start <- unname(interval[,"start"])
        
        for (m in seq_len(numProcesses)) {
          jobName <- paste0(prefix, "_", jobIter)
          jobIter <- jobIter + 1L
          
          end <- start + numReps[m] - 1L
                
          args <- paste0("'s/method/", methods[i], "/;s/overlap/", if (overlap[j]) "true" else "false",
                         "/;s/covariates/", covariates[k], "/;s/verbose/", if (verbose) "true" else "false",
                         "/;s/jobname/", jobName, "/;s/start/", start, "/;s/end/", end, "/'")
          system2("sed", args,
                  stdout = paste0("jobs/", jobName, ".q"),
                  stdin  = "template.q")
          system2("qsub", paste0("jobs/", jobName, ".q"))
          
          start <- end + 1L
        }
      }
    }
  }
}
