install_all_packages <- function(
        repo = "broadinstitute/bican-mccarroll-manuscript1",
        ref = "main",
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
        install.packages("remotes", repos = "https://cloud.r-project.org")
    }

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
    }

    invisible(NULL)
}

#install_all_packages()
#install_all_packages(ref = "jn_differential_expression")
#Rscript -e 'source("tools/install_packages.R"); install_all_packages(ref="jn_differential_expression")'

