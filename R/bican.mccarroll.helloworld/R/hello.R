# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   Document Package:          'Cmd + Shift + D'

#' Hello, world!
#' @param mom If TRUE, say hello to mom instead of world.
#' @export
hello <- function(mom=FALSE) {
  if (!is.logical(mom)) {
    stop("Argument 'mom' must be TRUE or FALSE.")
  } else if (mom) {
    print("Hi, mom!")
  } else {
    print("Hello, world!")
  }
}
