#' 	Print a summary of an output of the ecolRxC function
#'
#' @description Print method for objects obtained with the ecolRxC function.
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @param x An object output of the **ecolRxC** function.
#' @param ... Other arguments passed on to methods. Not currently used.
#' @param margins A TRUE/FALSE argument informing if the margins of the transition matrix should be displayed. Default TRUE.
#' @param digits Integer indicating the number of decimal places to be shown. Default, 2.
#'
#' @return
#' {No return value, called for side effects.} 
#'
#' @export
#' @method print ecolRxC
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
#' print(example, digits = 1, margins = TRUE)
#'

#' @export
print.ecolRxC  <- function(x,
                         ...,
                         margins = TRUE,
                         digits = 2)
{

  print.summary.ecolRxC(x = summary.ecolRxC(object = x),
                        margins = margins,
                        digits = digits)

}

#' 	Print a summary of a summary.ecolRxC object
#'
#' @description Print method for `summary.ecolRxC` objects
#' @inheritParams print.ecolRxC
#' @param x An `summary.ecolRxC` class object.
#' @return 
#' {No return value, called for side effects.} 
#' @method print summary.ecolRxC
#' @export
print.summary.ecolRxC  <- function(x,
                                 ...,
                                 margins = TRUE,
                                 digits = 2)
{

    tabla <- format(round(x$prop.matrix*100, digits), nsmall = digits)
    tabla <- apply(tabla, 2, as.character)
    rownames(tabla) <- rownames(x$prop.matrix)

    if (margins){
      nr <- nrow(tabla)
      tabla <- rbind(tabla, x$col.margins[1L:ncol(tabla)])
      tabla <- cbind(tabla,  c(format(x$row.margins[1L:nr], justify = "right"), ""))
    }

    cat("Estimated row-standardized transfer matrix \n")
    print(as.table(tabla))

}
