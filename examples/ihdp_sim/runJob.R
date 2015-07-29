start      <- as.integer(Sys.getenv("START"))
end        <- as.integer(Sys.getenv("END"))
overlap    <- if ((overlap <- Sys.getenv("OVERLAP")) == "") TRUE else as.logical(overlap)
covariates <- if ((covariates <- Sys.getenv("COVARIATES")) == "") "full" else covariates
method     <- Sys.getenv("METHOD")
verbose    <- if ((verbose <- Sys.getenv("VERBOSE")) == "") TRUE else as.logical(verbose)

require(npci)
source("results.R")

prefix <- getPrefix(method, overlap, covariates)

resultsFileName <- paste0(prefix, "_", start, "_", end, ".RData")
resultsFile <- file.path("data", resultsFileName)

if (file.exists(resultsFileName)) q("no")

numReps <- end - start + 1L

results <- matrix(NA, numReps, 6L)

colnames(results) <- 
  c("bias", "cov", "cil", "wrong", "tau.est", "precision")

source("data.R")

x <- as.matrix(x)
w <- rep(0.5, ncol(x))

prob.z <- glm(z ~ x, family = binomial)$fitted

estimand <- if (identical(overlap, TRUE)) "att" else "atc"

for (i in seq_len(numReps)){ 
  iter <- i + start - 1L
  
  if (verbose) cat("  running iter ", iter, "\n", sep = "")
  
  ## places mu.0, mu.1, y.0, y.1, and y into calling env
  generateDataForIterInCurrentEnvironment(iter, x, z, w, overlap, covariates)
    
  meanEffects <- if (identical(overlap, TRUE)) mu.1[z == 1] - mu.0[z == 1] else mu.1[z == 0] - mu.0[z == 0]
  results[i, "tau.est"] <- mean(y.1[z == 0] -  y.0[z == 0])
  
  if (method == "bart") {
    ## bart section
    
    ## each sample is a treatment effect, by person and by replication
    treatmentEffectSamples <- ci.estimate(y, x.r, z, method = "bart", estimand = estimand)
    ## average over people first, produces samples of ATC or ATT
    estimandSamples <- apply(treatmentEffectSamples, 2L, mean)
    results[i, "bias"] <- 4 - mean(estimandSamples)
    #ci <- quantile(atcSamples, c(0.025, 0.975))
    ci <- mean(estimandSamples) + sd(estimandSamples) * qnorm(c(0.025, 0.975))
    results[i, "cov"] <- if (ci[1L] < 4 && ci[2L] > 4) 1 else 0
    results[i, "cil"] <- ci[2L] - ci[1L]
    results[i, "wrong"] <- if (ci[1L] < 0 && ci[2L] > 0) 1 else 0
  
    ## average over samples, produces per-person TE estimates
    treatmentEffects <- apply(treatmentEffectSamples, 1L, mean)
    results[i, "precision"] <- sqrt(mean((treatmentEffects - meanEffects)^2))
    
  } else if (method == "bartipw") {
    treatmentEffectSamples <- ci.estimate(y, x.r, z, method = "bart", estimand = estimand, prob.z = prob.z)
    atcSamples <- apply(treatmentEffectSamples, 2L, mean)
    results[i, "bias"] <- 4 - mean(atcSamples)
    #ci <- quantile(atcSamples, c(0.025, 0.975))
    ci <- mean(atcSamples) + sd(atcSamples) * qnorm(c(0.025, 0.975))
    results[i, "cov"] <- if (ci[1L] < 4 && ci[2L] > 4) 1 else 0
    results[i, "cil"] <- ci[2L] - ci[1L]
    results[i, "wrong"] <- if (ci[1L] < 0 && ci[2L] > 0) 1 else 0
  
    treatmentEffects <- apply(treatmentEffectSamples, 1, mean)
    results[i, "precision"] <- sqrt(mean((treatmentEffects - meanEffects)^2))
  } else if (method %in% c("naive1", "naive2")) {
    atc <- ci.estimate(y, x.r, z, method = method, estimand = estimand)
    est <- mean(atc$te)
    se  <- atc$se
    results[i, "bias"] <- 4 - est
    ci <- est + se * qnorm(c(0.025, 0.975))
    results[i, "cov"] <- if (ci[1L] < 4 && ci[2L] > 4) 1 else 0
    results[i, "cil"] <- ci[2L] - ci[1L]
    results[i, "wrong"] <- if (ci[1L] < 0 && ci[2L] > 0) 1 else 0
    results[i, "precision"] <- sqrt(mean((atc$te - meanEffects)^2))
  } else {
    stop("method '", method, "' unrecognized")
  }
}

save(results, file = resultsFile)
