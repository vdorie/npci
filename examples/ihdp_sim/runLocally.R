numReps <- 250L

## this needs to be modified to adapt to work with results files already present
## from cluster runs

methods <- c("bart", "naive1", "naive2")

for (method in methods) {
  resultsFileName <- paste0(method, "_1_", numReps, ".RData")
  resultsFile <- file.path("data", resultsFileName)
  
  if (file.exists(resultsFile)) next
  
  Sys.setenv(START = 1)
  Sys.setenv(END = numReps)
  Sys.setenv(METHOD = method)
    
  source("runJob.R", local = TRUE, keep.source = FALSE)
}
