numReps <- 250L

methods <- "grouped"
#methods <- "bart"
#methods <- "btgp"
overlap <- c(TRUE, FALSE)
#covariates <- c(50, "junk", "select", "reduced")
covariates <- "reduced"
settings <- c("A", "B", "C")
p.score <- c("none", "true", "logistic", "bayesglm")
verbose <- TRUE

# covariates <- "full"
# p.score <- c("none", "true", "logistic")

baseResultInterval <- matrix(c(1L, numReps), 1L, 2L,
                             dimnames = list(NULL, c("start", "end")))

source("results.R")

for (i in seq_along(methods)) {
  for (j in seq_along(overlap)) {
    for (k in seq_along(covariates)) {
      for (l in seq_along(settings)) {
        for (m in seq_along(p.score)) {
          existingResults <- getResultIntervals(methods[i], overlap[j], covariates[k], settings[l], p.score[m])
          
          resultsIntervals <- intervalSubtraction(baseResultInterval, existingResults)
          if (nrow(resultsIntervals) == 0) next
          
          Sys.setenv(METHOD = methods[i])
          Sys.setenv(OVERLAP = if (overlap[j]) "true" else "false")
          Sys.setenv(COVARIATES = covariates[k])
          Sys.setenv(SETTING = settings[l])
          Sys.setenv(PROPENSITY = p.score[m])
          Sys.setenv(VERBOSE = verbose)
          
          for (n in seq_len(nrow(resultsIntervals))) {
            interval <- resultsIntervals[n,,drop = FALSE]
            #if (interval[,"start"] == 134L && covariates[k] == "reduced" && settings[l] == "B" && p.score[m] == "bayesglm") next
            Sys.setenv(START = interval[,"start"])
            Sys.setenv(END   = interval[,"end"])
            
            cat("running method '", methods[i], "' with overlap = ", overlap[j], ", covariates = '", covariates[k],
                "', setting = ", settings[l], ", prop = ", p.score[m], " from ", interval[,"start"], " to ", interval[,"end"], "\n", sep = "")
            
            source("runJob.R", local = new.env(parent = .GlobalEnv), keep.source = FALSE)
          }
        }
      }
    }
  }
}
