setwd("~/Repositories/npci/examples/ihdp_sim")

require(npci)

testIter <- 166L

source("results.R")

results.n1 <- (function() { res <- collateResults("naive1"); res$results[testIter,] })()
results.n2 <- (function() { res <- collateResults("naive2"); res$results[testIter,] })()
results.bt <- (function() { res <- collateResults("bart"); res$results[testIter,] })()

source("data.R")

x <- as.matrix(x)
w <- rep(0.5, ncol(x))

generateDataForIterInCurrentEnvironment(166L, x, z, w)

fitsFile <- file.path("data", "compFits.RData")
if (file.exists(fitsFile)) {
  load(fitsFile)
} else {
  fit.bt <- ci.fit(y, x, z, method = "bart")
  fit.n1 <- ci.fit(y, x, z, method = "naive1")
  fit.n2 <- ci.fit(y, x, z, method = "naive2")
  save(fit.bt, fit.n1, fit.n2, file = fitsFile)
}
rm(fitsFile)



colAverages <- sapply(seq_len(ncol(x)), function(j) {
  col <- x[,j]
  elm <- unique(col)
  if (length(elm) == 2) {
    temp <- col == elm[1]
    n1 <- sum(temp); n2 <- sum(!temp)
    return(if (n1 >= n2) elm[1] else elm[2])
  }
  mean(col)
})
names(colAverages) <- colnames(x)

columnIsContinuous <- apply(x, 2, function(col) is.numeric(col) && length(unique(col)) > 2)


numPoints <- 101
testColInd <- 5

xRange <- range(x[,testColInd])
xRange <- 1.02 * (xRange - mean(xRange)) + mean(xRange)
xValues <- seq(xRange[1], xRange[2], length.out = numPoints)


## 1st member - plot at average
## every subsequent member, go -0.5sd to +0.5sd for continuous
## and flip 0/1 for binary
lines.bt.0 <- matrix(NA, numPoints, 1 + 2 * (sum(columnIsContinuous) - 1) + sum(!columnIsContinuous))
lines.bt.1 <- lines.bt.0
lines.n1.0 <- lines.bt.0
lines.n1.1 <- lines.bt.0
lines.n2.0 <- lines.bt.0
lines.n2.1 <- lines.bt.0



x.test <- matrix(colAverages, numPoints, length(colAverages), byrow = TRUE)
x.test[,testColInd] <- xValues

pred.bt.0 <- predict(fit.bt, x.test, 0, keeptrainfits = FALSE)
pred.bt.1 <- predict(fit.bt, x.test, 1, keeptrainfits = FALSE)
pred.n1.0 <- predict(fit.n1, x.test, 0, "se")
pred.n1.1 <- predict(fit.n1, x.test, 1, "se")
pred.n2.0 <- predict(fit.n2, x.test, 0, "se")
pred.n2.1 <- predict(fit.n2, x.test, 1, "se")


lines.bt.0[,1] <- pred.bt.0$mean
lines.bt.1[,1] <- pred.bt.1$mean
lines.n1.0[,1] <- pred.n1.0$fit
lines.n1.1[,1] <- pred.n1.1$fit
lines.n2.0[,1] <- pred.n2.0$fit
lines.n2.1[,1] <- pred.n2.1$fit

poly.bt.0 <- cbind(pred.bt.0$bound[1,], rev(pred.bt.0$bound[2,]))
poly.bt.1 <- cbind(pred.bt.1$bound[1,], rev(pred.bt.1$bound[2,]))
poly.n1.0 <- cbind(pred.n1.0$fit - pred.n1.0$se * qnorm(0.975), rev(pred.n1.0$fit + pred.n1.0$se * qnorm(0.975)))
poly.n1.1 <- cbind(pred.n1.1$fit - pred.n1.1$se * qnorm(0.975), rev(pred.n1.1$fit + pred.n1.1$se * qnorm(0.975)))
poly.n2.0 <- cbind(pred.n2.0$fit - pred.n2.0$se * qnorm(0.975), rev(pred.n2.0$fit + pred.n2.0$se * qnorm(0.975)))
poly.n2.1 <- cbind(pred.n2.1$fit - pred.n2.1$se * qnorm(0.975), rev(pred.n2.1$fit + pred.n2.1$se * qnorm(0.975)))


lineIndex <- 2
for (colInd in seq_len(ncol(x))) {
  if (colInd == testColInd) next
  
  if (columnIsContinuous[colInd]) 
    x.test[,colInd] <- colAverages[colInd] - 0.5 * sd(x[,colInd])
  else
    x.test[,colInd] <- 1 - x.test[1,colInd]
  
  lines.bt.0[,lineIndex] <- predict(fit.bt, x.test, 0, summarize = quote(mean(x)), keeptrainfits = FALSE)
  lines.bt.1[,lineIndex] <- predict(fit.bt, x.test, 1, summarize = quote(mean(x)), keeptrainfits = FALSE)
  lines.n1.0[,lineIndex] <- predict(fit.n1, x.test, 0)
  lines.n1.1[,lineIndex] <- predict(fit.n1, x.test, 1)
  lines.n2.0[,lineIndex] <- predict(fit.n2, x.test, 0, error = "se")$fit
  lines.n2.1[,lineIndex] <- predict(fit.n2, x.test, 1, error = "se")$fit
  
  
  lineIndex <- lineIndex + 1
  
  
  if (columnIsContinuous[colInd]) {
    x.test[,colInd] <- colAverages[colInd] + 0.5 * sd(x[,colInd])
    
    lines.bt.0[,lineIndex] <- predict(fit.bt, x.test, 0, summarize = quote(mean(x)), keeptrainfits = FALSE)
    lines.bt.1[,lineIndex] <- predict(fit.bt, x.test, 1, summarize = quote(mean(x)), keeptrainfits = FALSE)
    lines.n1.0[,lineIndex] <- predict(fit.n1, x.test, 0)
    lines.n1.1[,lineIndex] <- predict(fit.n1, x.test, 1)
    lines.n2.0[,lineIndex] <- predict(fit.n2, x.test, 0, error = "se")$fit
    lines.n2.1[,lineIndex] <- predict(fit.n2, x.test, 1, error = "se")$fit
    
    lineIndex <- lineIndex + 1
  }
  
  x.test[,colInd] <- colAverages[colInd]
}


wrapper.0 <- function(x) {
  x.test <- c(1, x.test)
  sapply(x, function(x.i) {
    x.test[colInd] <- x.i
    f.0(x.test)[1]
  })
}
wrapper.1 <- function(x) {
  x.test <- c(1, x.test)
  sapply(x, function(x.i) {
    x.test[colInd] <- x.i
    f.1(x.test)[1]
  })
}
environment(wrapper.0) <- npci:::args2env(baseenv(), colInd = testColInd + 1, f.0, f.1, x.test = colAverages)
environment(wrapper.1) <- environment(wrapper.0)

yRange <- range(y,
                lines.bt.0, lines.bt.1, poly.bt.0, poly.bt.1,
                lines.n1.0, lines.n1.1, poly.n1.0, poly.n1.1,
                lines.n2.0, lines.n2.1, poly.n2.0, poly.n2.1)
yRange <- 1.02 * (yRange - mean(yRange)) + mean(yRange)

poly.x <- cbind(xValues, rev(xValues))




pdf("~/Desktop/gp_comp.pdf", width = 6, height = 2.5)
par(mfrow = c(1, 3))
plot(NULL, type = "n", xlim = xRange, ylim = yRange, xlab = colnames(x)[testColInd],
     ylab = "Response", main = "BART")

polygon(poly.x, poly.bt.0, col = rgb(0.95, 0.95, 0.95), border = NA)
polygon(poly.x, poly.bt.1, col = rgb(0.90, 0.90, 0.90), border = NA)

for (lineIndex in seq(2L, ncol(lines.bt.0))) {
  lines(xValues, lines.bt.0[,lineIndex], col = "gray", lwd = 0.5)
  lines(xValues, lines.bt.1[,lineIndex], col = "black", lwd = 0.5)
}
points(x[,testColInd], y, col = ifelse(z == 1, "black", "gray"), pch = 20, cex = 0.6)

lines(xValues, lines.bt.0[,1], col = "blue", lwd = 1)
lines(xValues, lines.bt.1[,1], col = "blue", lwd = 1)

curve(wrapper.0, add = TRUE, col = "red")
curve(wrapper.1, add = TRUE, col = "red")


plot(NULL, type = "n", xlim = xRange, ylim = yRange, xlab = colnames(x)[testColInd],
     ylab = "Response", main = "Naive 1")

polygon(poly.x, poly.n1.0, col = rgb(0.95, 0.95, 0.95), border = NA)
polygon(poly.x, poly.n1.1, col = rgb(0.90, 0.90, 0.90), border = NA)

for (lineIndex in seq(2L, ncol(lines.n1.0))) {
  lines(xValues, lines.n1.0[,lineIndex], col = "gray", lwd = 0.5)
  lines(xValues, lines.n1.1[,lineIndex], col = "black", lwd = 0.5)
}
points(x[,testColInd], y, col = ifelse(z == 1, "black", "gray"), pch = 20, cex = 0.6)

lines(xValues, lines.n1.0[,1], col = "blue", lwd = 1)
lines(xValues, lines.n1.1[,1], col = "blue", lwd = 1)

curve(wrapper.0, add = TRUE, col = "red")
curve(wrapper.1, add = TRUE, col = "red")


plot(NULL, type = "n", xlim = xRange, ylim = yRange, xlab = colnames(x)[testColInd],
     ylab = "Response", main = "Naive 2")

polygon(poly.x, poly.n2.0, col = rgb(0.95, 0.95, 0.95), border = NA)
polygon(poly.x, poly.n2.1, col = rgb(0.90, 0.90, 0.90), border = NA)

for (lineIndex in seq(2L, ncol(lines.n2.0))) {
  lines(xValues, lines.n2.0[,lineIndex], col = "gray", lwd = 0.5)
  lines(xValues, lines.n2.1[,lineIndex], col = "black", lwd = 0.5)
}
points(x[,testColInd], y, col = ifelse(z == 1, "black", "gray"), pch = 20, cex = 0.6)

lines(xValues, lines.n2.0[,1], col = "blue", lwd = 1)
lines(xValues, lines.n2.1[,1], col = "blue", lwd = 1)

curve(wrapper.0, add = TRUE, col = "red")
curve(wrapper.1, add = TRUE, col = "red")

dev.off()
