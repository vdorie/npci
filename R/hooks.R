.onUnload <- function(libpath)
{
  if (is.loaded("npci_grouped_updateCovMatrix", PACKAGE = "npci")) {
    library.dynam.unload("npci", libpath)
  }
}

