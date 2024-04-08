.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s. This is due to a conflict with registered trademark of company.",
    pkgname, "3.19")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}

