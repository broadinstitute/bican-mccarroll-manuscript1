install_all_packages <- function(
        repo = "broadinstitute/bican-mccarroll-manuscript1",
        ref = "main",
        path = NULL,
        pkgs = c(
            "R/bican.mccarroll.differentialexpression",
            "R/bican.mccarroll.de.analysis",
            "R/bican.mccarroll.eqtl",
            "R/bican.mccarroll.figures"
        ),
        dependencies = TRUE,
        upgrade = "never",
        force = FALSE) {

    if (!requireNamespace("remotes", quietly = TRUE)) {
        # Use https://cloud.r-project.org to install remotes only if the default CRAN mirror is not set.
        if (getOption("repos")["CRAN"] == "@CRAN@") {
            install.packages("remotes", repos = "https://cloud.r-project.org")
        } else {
            install.packages("remotes")
        }
    }

    if (is.null(path)) {
        for (subdir in pkgs) {
            message("Installing ", subdir, " from branch/ref: ", ref)
            remotes::install_github(
                repo = repo,
                subdir = subdir,
                ref = ref,
                dependencies = dependencies,
                upgrade = upgrade,
                force = force
            )
            # exit now if the library failed to install
            library(basename(subdir), character.only = TRUE)
        }
    } else {
        for (subdir in pkgs) {
            message("Installing ", subdir, " from local path: ", path)
            remotes::install_local(
                path = file.path(path, subdir),
                dependencies = dependencies,
                upgrade = upgrade,
                force = force
            )
            # exit now if the library failed to install
            library(basename(subdir), character.only = TRUE)
        }
    }

    invisible(NULL)
}

#install_all_packages()
#install_all_packages(ref = "jn_differential_expression")
#Rscript -e 'source("tools/install_packages.R"); install_all_packages(ref="jn_differential_expression")'
