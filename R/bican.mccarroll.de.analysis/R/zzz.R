# Data.table "aware" package loading
.onLoad <- function(libname, pkgname) {
    invisible(NULL)
}
utils::globalVariables(":=")
