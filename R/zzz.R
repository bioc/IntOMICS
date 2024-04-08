.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed and excluded from  Bioconductor version %s due to conflicting trademark licensing.",
    pkgname, "3.19")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}

