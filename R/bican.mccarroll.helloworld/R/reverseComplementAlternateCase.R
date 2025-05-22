
#' Reverse Complement of a DNA Sequence with Alternating Case
#' E.g.: "ATCG" -> "CgAt"
#' Demonstrates the used of a dependency.
#' @param s A string representing a DNA sequence.
#' @return A string representing the reverse complement of the input DNA sequence, with alternating case.
#' @export
#' @import spgs
reverseComplementAlternateCase<-function(s) {
  rc = spgs::reverseComplement(s)
  chars = strsplit(rc, "")[[1]]
  for (i in seq_along(chars)) {
    if (i %% 2 == 1) {
      chars[i] = toupper(chars[i])
    } else {
      chars[i] = tolower(chars[i])
    }
  }
  return(paste(chars, collapse = ""))
}
