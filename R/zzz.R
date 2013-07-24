.packageName <- "inapplicable"

.onLoad <- function(libname, pkgname) {
    library.dynam("fitch", pkgname, libname)
}

