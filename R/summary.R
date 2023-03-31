#' 	Summarize a ecolRxC output object
#'
#' @description Summary method for objects obtained with the ecolRxC function
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @param object An object output of the **ecolRxC** function.
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @return
#' An object of class `"summary.ecolRxC"`. 
#' A list with four components:
#'  \item{prop.matrix}{ A matrix of order JxK (RxC) with the estimated proportions of the row-standardized vote transitions from election 1 to election 2.}
#'  \item{counts.matrix}{ A matrix of order JxK (RxC) with the estimated vote transfers from election 1 to election 2.}
#'  \item{row.margins}{ A vector of length R with aggregate observed distribution of votes in election 1.}
#'  \item{col.margins}{ A vector of length C with aggregate observed distribution of votes in election 2.}
#'
#' @export
#' @method summary ecolRxC
#' @examples
#' votes1 <- structure(list(P1 = c(16L, 4L, 13L, 6L, 1L, 16L, 6L, 17L, 48L, 14L),
#'                          P2 = c(8L, 3L, 0L, 5L, 1L, 4L, 7L, 6L, 28L, 8L),
#'                          P3 = c(38L, 11L, 11L, 3L, 13L, 39L, 14L, 34L, 280L, 84L),
#'                          P4 = c(66L, 5L, 18L, 39L, 30L, 57L, 35L, 65L, 180L, 78L),
#'                          P5 = c(14L, 0L, 5L, 2L, 4L, 21L, 6L, 11L, 54L, 9L),
#'                          P6 = c(8L, 2L, 5L, 3L, 0L, 7L, 7L, 11L, 45L, 17L),
#'                          P7 = c(7L, 3L, 5L, 2L, 3L, 17L, 7L, 13L, 40L, 8L)),
#'                          row.names = c(NA, 10L), class = "data.frame")
#' votes2 <- structure(list(C1 = c(2L, 1L, 2L, 2L, 0L, 4L, 0L, 4L, 19L, 14L),
#'                          C2 = c(7L, 3L, 1L, 7L, 2L, 5L, 3L, 10L, 21L, 6L),
#'                          C3 = c(78L, 7L, 28L, 42L, 28L, 84L, 49L, 85L, 260L, 100L),
#'                          C4 = c(56L, 14L, 20L, 7L, 19L, 54L, 22L, 50L, 330L, 91L),
#'                          C5 = c(14L, 3L, 6L, 2L, 3L, 14L, 8L, 8L, 45L, 7L)),
#'                          row.names = c(NA, 10L), class = "data.frame")
#' example <- ecolRxC(votes1, votes2, method = "IPF")
#' summary(example)
#'

summary.ecolRxC <- function(object, ...){

  output <- NULL
  output$prop.matrix <- object$VTM
  output$counts.matrix <- object$VTM.votes

  output$row.margins <- round(rowSums(object$VTM.votes))
  output$col.margins <- round(colSums(object$VTM.votes))

  class(output) <- "summary.ecolRxC"
  output

}
