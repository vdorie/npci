## WARNING

**This package is experimental and not supported. To fit a non-parametric model for causal inference, use [bartCause](https://github.com/vdorie/bartCause) intead.**

## Prequisites

[Rtools](https://cran.r-project.org/bin/windows/Rtools/) for windows. [gfortran](https://github.com/fxcoudert/gfortran-for-macOS/releases) for Mac OS.

## Notes

Contains R code to fit non-parmetric causal inference models. Of particular interest are the calculations in the examples folder.

To run the examples, first install the package by running, from within R

    install.packages("path/to/repository", type = "source", repos = NULL)

or from the command line

    R CMD INSTALL path/to/repository

The file `examples/toy/generatePlots.R` contains code to produce the toy data example. The simulation data can be simply loaded by running `examples/ihdp_sim/loadData.R`, or the simulations produced by running either `examples/ihdp_sim/runLocally.R` or `examples/ihdp_sim/queueJobs.R` (if on a cluster using TORQUE). Change your working directory before doing any of the above.
