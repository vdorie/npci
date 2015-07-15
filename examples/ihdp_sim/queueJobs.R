## for use with TORQUE/PBS/qsub

totalNumReps <- 250L

methods <- c("bart", "naive1", "naive2")
numRepsPerProcess <- c(250L, 5L, 10L)
setting <- "lowp"

if (!dir.exists("jobs")) dir.create("jobs")

baseResultInterval <- matrix(c(1L, totalNumReps), 1L, 2L,
                             dimnames = list(NULL, c("start", "end")))

source("results.R")

jobIter <- 1L

for (i in seq_along(methods)) {
  method <- methods[i]
  prefix <- if (is.null(setting)) method else paste0(method, "_", setting)
  
  existingResults <- getResultIntervals(method, setting)
  
  resultsIntervals <- intervalSubtraction(baseResultInterval, existingResults)
  if (nrow(resultsIntervals) == 0) next
  
  for (j in seq_len(nrow(resultsIntervals))) {
    interval <- resultsIntervals[j,, drop = FALSE]
    numRepsPerInterval <- unname(interval[,"end"] - interval[,"start"] + 1L)
    
    numProcesses <- numRepsPerInterval %/% numRepsPerProcess[i] +
      if (numRepsPerInterval %% numRepsPerProcess[i] != 0) 1L else 0L
    numReps <- rep(numRepsPerInterval %/% numProcesses, numProcesses) +
      c(rep(1L, numRepsPerInterval %% numProcesses), rep(0L, numProcesses - numRepsPerInterval %% numProcesses))
    
    start <- unname(interval[,"start"])
    
    for (k in seq_len(numProcesses)) {
      jobName <- paste0(prefix, "_", jobIter)
      jobIter <- jobIter + 1L
      
      end <- start + numReps[k] - 1L
            
      args <- paste0("'s/jobname/", jobName, "/;s/method/", method, "/;s/setting/", setting, "/;s/start/", start, "/;s/end/", end, "/'")
      system2("sed", args,
              stdout = paste0("jobs/", jobName, ".q"),
              stdin  = "template.q")
      system2("qsub", paste0("jobs/", jobName, ".q"))
      
      start <- end + 1L
    }
  }
}
