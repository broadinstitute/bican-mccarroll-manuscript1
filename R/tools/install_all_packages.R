install_all_packages <- function(
        repo = "broadinstitute/bican-mccarroll-manuscript1",
        ref = "main",
        pkgs = c(
            "R/bican.mccarroll.differentialexpression",
            "R/bican.mccarroll.eqtl",
            "R/bican.mccarroll.figures"
        ),
        dependencies = TRUE,
        upgrade = "never") {

    if (!requireNamespace("remotes", quietly = TRUE)) {
        stop("Please install remotes: install.packages('remotes')")
    }

    for (subdir in pkgs) {
        message("Installing ", subdir, " from branch/ref: ", ref)
        remotes::install_github(
            repo = repo,
            subdir = subdir,
            ref = ref,
            dependencies = dependencies,
            upgrade = upgrade
        )
    }

    invisible(NULL)
}

install_all_packages()

