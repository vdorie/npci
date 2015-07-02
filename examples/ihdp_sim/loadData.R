dataFile <- file.path("data", "ihdp.RData")
if (!file.exists(dataFile)) stop("ihdp data file not available")

load(dataFile)

ihdp <- subset(ihdp, treat != 1 | momwhite != 0)

covariateNames <- c("bw", "b.head", "preterm", "birth.o", "nnhealth", "momage",
                    "sex", "twin", "b.marr", "mom.lths", "mom.hs", "mom.scoll",
                    "cig", "first", "booze", "drugs", "work.dur", "prenatal",
                    "ark", "ein", "har", "mia", "pen", "tex", "was")

x <- ihdp[, covariateNames]
trans <- npci:::getTransformations(x)

x <- npci:::transform(x, trans$standardize)
z <- ihdp$treat

rm(dataFile, ihdp, covariateNames, trans)
