install_all_packages <- function(repo = "broadinstitute/bican-mccarroll-manuscript1",
                                      ref = "main",
                                      pkgs = c(
                                          "R/bican.mccarroll.differentialexpression",
                                          "R/bican.mccarroll.eqtl",
                                          "R/bican.mccarroll.figures"
                                      ),
                                      upgrade = "never",
                                      dependencies = TRUE) {

    if (!requireNamespace("remotes", quietly = TRUE)) {
        stop("Please install remotes: install.packages('remotes')")
    }

    for (subdir in pkgs) {
        remotes::install_github(
            repo,
            ref = ref,
            subdir = subdir,
            dependencies = dependencies,
            upgrade = upgrade
        )
    }

    invisible(NULL)
}

install_all_packages()

