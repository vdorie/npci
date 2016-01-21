## for use with TORQUE/PBS/qsub

if (require(npci, quietly = TRUE) == FALSE)
  stop("npci package not available")

totalNumReps <- 250L

#methods <- c("grouped", "independent", "bart")
methods <- "grouped"
numRepsPerProcess <- 5L
overlap <- c(TRUE, FALSE)
covariates <- c("full", 50, "select", "reduced", "junk")
settings <- c("A", "B", "C")
p.score <- c(TRUE, FALSE)
verbose <- TRUE

if (!dir.exists("jobs")) dir.create("jobs")

baseResultInterval <- matrix(c(1L, totalNumReps), 1L, 2L,
                             dimnames = list(NULL, c("start", "end")))

source("results.R")

jobIter <- 1L

for (i in seq_along(methods)) {
  for (j in seq_along(overlap)) {
    for (k in seq_along(covariates)) {
      for (l in seq_along(settings)) {
        for (m in seq_along(p.score)) {
          
          prefix <- getPrefix(methods[i], overlap[j], covariates[k], settings[l], p.score[m])
  
          existingResults <- getResultIntervals(methods[i], overlap[j], covariates[k], settings[l], p.score[m])
          
          resultsIntervals <- intervalSubtraction(baseResultInterval, existingResults)
          if (nrow(resultsIntervals) == 0) next
          
          for (n in seq_len(nrow(resultsIntervals))) {
            interval <- resultsIntervals[n,,drop = FALSE]
            numRepsPerInterval <- unname(interval[,"end"] - interval[,"start"] + 1L)
            
            numProcesses <- numRepsPerInterval %/% numRepsPerProcess[i] +
              if (numRepsPerInterval %% numRepsPerProcess[i] != 0) 1L else 0L
            numReps <- rep(numRepsPerInterval %/% numProcesses, numProcesses) +
              c(rep(1L, numRepsPerInterval %% numProcesses), rep(0L, numProcesses - numRepsPerInterval %% numProcesses))
            
            start <- unname(interval[,"start"])
            
            for (o in seq_len(numProcesses)) {
              jobName <- paste0(prefix, "_", jobIter)
              jobIter <- jobIter + 1L
              
              end <- start + numReps[o] - 1L
              
              args <- paste0("'s/_METHOD_/", methods[i], "/;s/_OVERLAP_/", if (overlap[j]) "true" else "false",
                             "/;s/_COVARIATES_/", covariates[k], "/;s/_SETTING_/", settings[l],
                             "/;s/_PROPENSITY_/", p.score[m], "/;s/_VERBOSE_/", if (verbose) "true" else "false",
                             "/;s/_JOBNAME_/", jobName, "/;s/_START_/", start, "/;s/_END_/", end, "/'")
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
  }
}
