#' Font used for all manuscript figures.
MANUSCRIPT_FONT <- "Arial"

#' Open an SVG device that renders text in the manuscript font
#'
#' Wraps \code{svglite::svglite()}, binding \code{MANUSCRIPT_FONT} to the
#' \code{"sans"} font alias so that ggplot2's default (unset) text family
#' resolves to it. Suitable both as a direct device opener
#' (\code{svglite_manuscript(file, width, height)}) and as a \code{ggsave}
#' \code{device} argument.
#'
#' @param filename Output SVG file path.
#' @param ... Additional arguments passed to \code{svglite::svglite()}.
#' @return Invisibly, the result of \code{svglite::svglite()}.
svglite_manuscript <- function(filename, ...) {
  svglite::svglite(
    filename,
    ...,
    system_fonts = list(sans = MANUSCRIPT_FONT)
  )
}

#' Warn if the manuscript font is not resolvable on this system
#'
#' Called at package load. \code{svglite_manuscript()} silently falls back to
#' a substitute font if \code{MANUSCRIPT_FONT} is not installed, so this check
#' exists to surface that condition loudly instead of letting figures
#' regenerate in the wrong font unnoticed.
.check_manuscript_font <- function() {
  match <- systemfonts::match_fonts(MANUSCRIPT_FONT)
  if (!grepl(MANUSCRIPT_FONT, basename(match$path), ignore.case = TRUE)) {
    warning(
      sprintf(
        "Manuscript font '%s' not found on this system; figures will fall back to '%s'.",
        MANUSCRIPT_FONT,
        basename(match$path)
      ),
      call. = FALSE
    )
  }
}
