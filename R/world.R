#' Function for print other message
#'
#' @param a a message
#' @param b a message
#'
#' @return Print word \code{a} and \code{b}
#'
#' @examples
#' world('Hello', 'world')
#'
#' @export
world <- function(message1, message2) {
  return(paste(message1, message2, sep=" "))
}
