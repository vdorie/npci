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


plotFit <- function(y, x, z, f0, f1, method, estimand = NULL, prob.z = NULL, main = method)
{
  fit <- ci.fit(y, x, z, method, estimand, prob.z)
  
  xRange <- range(x)
  xRange <- 1.02 * (xRange - mean(xRange)) + mean(xRange)
  xValues <- seq(xRange[1], xRange[2], length.out = 101)
  
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
  }
  
  poly.x <- cbind(xValues, rev(xValues))
  
  yRange <- range(poly.y.0, poly.y.1)
  yRange <- 1.02 * (yRange - mean(yRange)) + mean(yRange)
  
  plot(NULL, type = "n", xlim = range(xValues), ylim = yRange, xlab = "", ylab = "",
       main = main)
  polygon(poly.x, poly.y.0, col = rgb(0.95, 0.95, 0.95), border = NA)
  polygon(poly.x, poly.y.1, col = rgb(0.9, 0.9, 0.9), border = NA)
  lines(xValues, y.hat.1)
  lines(xValues, y.hat.0, col = "gray")
  lines(xValues, f0(xValues), lwd = 0.7, col = "red")
  lines(xValues, f1(xValues), lwd = 0.7, col = "red")
  points(x, y, col = ifelse(z == 1, "black", "gray"), pch = 20)
}

generatePlots <- function(path = ".") {
  toyData <- generateToyData()
  attach(toyData)
  
  pdf(file.path(path, "toy_naive1.pdf"), 3.5, 3.5)
  plotFit(y, x, z, f0, f1, "naive1")
  dev.off()
  
  pdf(file.path(path, "toy_naive2.pdf"), 3.5, 3.5)
  plotFit(y, x, z, f0, f1, "naive2")
  dev.off()
  
  pdf(file.path(path, "toy_bart.pdf"), 3.5, 3.5)
  par(mfrow = c(2, 2))
  plotFit(y, x, z, f0, f1, "bart", main = "No Weights")
  prob.z <- glm(z ~ x, family = binomial())$fitted
  plotFit(y, x, z, f0, f1, "bart", "ate", prob.z, "ATE Weights")
  plotFit(y, x, z, f0, f1, "bart", "att", prob.z, "ATT Weights")
  plotFit(y, x, z, f0, f1, "bart", "atc", prob.z, "ATC Weights")
  dev.off()
  
  detach(toyData)
}

generatePlots()

rm(generateToyData, plotFit, generatePlots)
