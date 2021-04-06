require(npci)

generateToyData <- function() {
  f0 <- function(x) 72 + 3 * sqrt(x)
  f1 <- function(x) 90 + exp(0.06 * x)
  
  oldSeed <- if (exists(".Random.seed")) .Random.seed else NULL
  set.seed(997L)
  
  ## generate true values
  n <- 120L
  p <- 0.5
  
  z <- rbinom(n, 1L, p)
  n1 <- sum(z); n0 <- n - n1
  
  x <- numeric(n)
  x[z == 0] <- rnorm(n0, 20, 10)
  x[z == 1] <- rnorm(n1, 40, 10)
  y <- numeric(n)
  y[z == 0] <- f0(x[z == 0])
  y[z == 1] <- f1(x[z == 1])
  y <- y + rnorm(n)
  
  if (!is.null(oldSeed)) .Random.seed <- oldSeed
  
  npci:::namedList(x, y, n, p, z, f0, f1)
}

getPlotXValues <- function(x, n.points = 101L)
{
  xRange <- range(x)
  xRange <- 1.02 * (xRange - mean(xRange)) + mean(xRange)
  seq(xRange[1], xRange[2], length.out = n.points)
}

getPlotYValues <- function(y, x, z, f0, f1, xValues, method, estimand = NULL, prob.z = NULL)
{
  if (method == "truth")
    return(list(y.0 = f0(xValues), y.1 = f1(xValues)))
  
  fit <- ci.fit(y, x, z, method, estimand, prob.z)
  
  if (method %in% c("naive1", "naive2")) {
    pred.0 <- predict(fit, xValues, 0, "se")
    pred.1 <- predict(fit, xValues, 1, "se")
    
    y.hat.0 <- pred.0$fit
    y.hat.1 <- pred.1$fit
    
    se.0    <- pred.0$se
    se.1    <- pred.1$se
    
    poly.y.0 <- cbind(y.hat.0 - se.0 * qnorm(0.975), rev(y.hat.0 + se.0 * qnorm(0.975)))
    poly.y.1 <- cbind(y.hat.1 - se.1 * qnorm(0.975), rev(y.hat.1 + se.1 * qnorm(0.975)))
  } else if (method == "bart") {
    pred.0 <- predict(fit, xValues, 0, keeptrainfits = FALSE)
    pred.1 <- predict(fit, xValues, 1, keeptrainfits = FALSE) 
    
    y.hat.0 <- pred.0$mean
    y.hat.1 <- pred.1$mean
    
    lower.0 <- pred.0$bound[1,]
    lower.1 <- pred.1$bound[1,]
    upper.0 <- pred.0$bound[2,]
    upper.1 <- pred.1$bound[2,]
    
    poly.y.0 <- cbind(upper.0, rev(lower.0))
    poly.y.1 <- cbind(upper.1, rev(lower.1))
  } else if (method %in%  c("grouped", "independent")) {
    optPars <- if (method == "grouped") fit@env$opt[[1]]$par else rbind(fit@env$opt.0[[1]]$par, fit@env$opt.1[[1]]$par)
    pred.0 <- predict(fit, xValues, 0, "se", pars = optPars)
    pred.1 <- predict(fit, xValues, 1, "se", pars = optPars)
    
    y.hat.0 <- pred.0$fit
    y.hat.1 <- pred.1$fit
    
    se.0    <- pred.0$se
    se.1    <- pred.1$se
    
    poly.y.0 <- cbind(y.hat.0 - se.0 * qnorm(0.975), rev(y.hat.0 + se.0 * qnorm(0.975)))
    poly.y.1 <- cbind(y.hat.1 - se.1 * qnorm(0.975), rev(y.hat.1 + se.1 * qnorm(0.975)))
  } else if (method == "btgp") {
    pred.0 <- predict(fit, xValues, 0, pred.n = FALSE)
    pred.1 <- predict(fit, xValues, 1, pred.n = FALSE)
    
    y.hat.0 <- pred.0$fit
    y.hat.1 <- pred.1$fit
    
    se.0    <- sqrt(diag(pred.0$vcov))
    se.1    <- sqrt(diag(pred.1$vcov))
    
    poly.y.0 <- cbind(y.hat.0 - se.0 * qnorm(0.975), rev(y.hat.0 + se.0 * qnorm(0.975)))
    poly.y.1 <- cbind(y.hat.1 - se.1 * qnorm(0.975), rev(y.hat.1 + se.1 * qnorm(0.975)))
  }
  
  list(y.0 = f0(xValues), y.hat.0 = y.hat.0, y.1 = f1(xValues), y.hat.1 = y.hat.1, poly.y.0 = poly.y.0, poly.y.1 = poly.y.1)
}

plotFit <- function(y, x, z, yValues, xValues, yRange, main)
{
  par(mar = c(1.5, 1.5, 1.4, 0.1))
  plot(NULL, type = "n", xlim = range(xValues), ylim = yRange, xlab = "", ylab = "",
       main = main)
  
  if (is.null(yValues$y.hat.0)) {
    lines(xValues, yValues$y.0, col = "gray")
    lines(xValues, yValues$y.1)
    points(x, y, col = ifelse(z == 1, "black", "gray"), pch = 20)
  } else {
    poly.x <- cbind(xValues, rev(xValues))
    
    polygon(poly.x, yValues$poly.y.0, col = rgb(0.95, 0.95, 0.95), border = NA)
    polygon(poly.x, yValues$poly.y.1, col = rgb(0.9, 0.9, 0.9), border = NA)
    lines(xValues, yValues$y.hat.0, col = "gray")
    lines(xValues, yValues$y.hat.1)
    lines(xValues, yValues$y.0, lwd = 0.7, col = "red")
    lines(xValues, yValues$y.1,, lwd = 0.7, col = "red")
    points(x, y, col = ifelse(z == 1, "black", "gray"), pch = 20)
  }
}

generatePlots <- function(path = ".") {
  toyData <- generateToyData()
  attach(toyData)
  
  xValues <- getPlotXValues(x)

  yValues.truth <- getPlotYValues(y, x, z, f0, f1, xValues, "truth")
  yValues.bart  <- getPlotYValues(y, x, z, f0, f1, xValues, "bart")
  yValues.group <- getPlotYValues(y, x, z, f0, f1, xValues, "grouped")
  yValues.indep <- getPlotYValues(y, x, z, f0, f1, xValues, "independent")
  yValues.btgp  <- getPlotYValues(y, x, z, f0, f1, xValues, "btgp")

  yRange <- range(yValues.truth$y.0, yValues.truth$y.1,
                  yValues.bart$y.hat.0, yValues.bart$y.hat.1,
                  yValues.group$y.hat.0, yValues.group$y.hat.1,
                  yValues.indep$y.hat.0, yValues.indep$y.hat.1,
                  yValues.btgp$y.hat.0, yValues.btgp$y.hat.1)
  yRange <- 1.02 * (yRange - mean(yRange)) + mean(yRange)

  
  pdf(file.path(path, "toy_truth.pdf"), 3.5, 3.5)
  plotFit(y, x, z, yValues.truth, xValues, yRange, "Truth")
  dev.off()
  
  pdf(file.path(path, "toy_bart.pdf"), 3.5, 3.5)
  plotFit(y, x, z, yValues.bart, xValues, yRange, "BART")
  dev.off()
  
  pdf(file.path(path, "toy_grouped.pdf"), 3.5, 3.5)
  plotFit(y, x, z, yValues.group, xValues, yRange, "GP-Grouped")
  dev.off()
  
  pdf(file.path(path, "toy_independent.pdf"), 3.5, 3.5)
  plotFit(y, x, z, yValues.indep, xValues, yRange, "GP-Independent")
  dev.off()
  
  pdf(file.path(path, "toy_btgp.pdf"), 3.5, 3.5)
  plotFit(y, x, z, yValues.btgp, xValues, yRange, "BTGP")
  dev.off()
  
  detach(toyData)
}

generatePlots()

rm(generateToyData, plotFit, generatePlots)
