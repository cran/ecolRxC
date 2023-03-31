#' Ecological Inference of RxC Tables by Latent Structure Approaches
#'
#' @description  Estimates JxK (RxC) vote transfer matrices (ecological contingency tables) based on Thomsen (1987) and Park (2008) approaches.
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @references Achen, C.H. (2000). The Thomsen Estimator for Ecological Inference (Unpublished manuscript). University of Michigan.
#' @references Park, W.-H. (2008). Ecological Inference and Aggregate Analysis of Elections. PhD Dissertation. University of Michigan.
#' @references Pavia, J.M. (2022). Adjustment of initial estimates of voter transition probabilities to guarantee consistency and completeness.
#' @references Thomsen, S.R. (1987). Danish Elections 1920-79: a Logit Approach to Ecological Analysis and Inference. Politica, Aarhus, Denmark.
#'
#' @param votes.election1 data.frame (or matrix) of order IxJ1 with the votes gained by
#'                        (or the counts corresponding to) the J1 (social classes) political options competing
#'                        (available) on election 1 (or origin) in the I units considered.
#'
#' @param votes.election2 data.frame (or matrix) of order IxK2 with the votes gained by
#'                        (or the counts corresponding to) the K2 political options competing
#'                        (available) on election 2 (or destination) in the I (territorial) units considered.
#'
#' @param scale A character string indicating the type of transformation to be applied to the vote
#'              fractions for applying ecological inference. Only `logit` and `probit` are allowed.
#'              Default, `probit`.
#'
#' @param method A character string indicating the algorithm to be used for adjusting (making congruent with
#'               the observed margins) the initial crude fractions attained in a 2x2 fashion.
#'               Only `Thomsen` (see sec. 4.3 in Thomsen, 1987) and `IPF` (iterative proportional fitting,
#'               also known as raking) are allowed. This argument has no effect in the 2x2 case.
#'               Default, `Thomsen`.
#'
#' @param local A TRUE/FALSE argument indicating whether local solutions (solutions for each polling unit)
#'              must be computed. In that case, the global solution is attained as composition/aggregation
#'              of local solutions. When `method = "Thomsen"` local solutions are always computed.
#'              Default `TRUE`.
#'
#' @param census.changes A character string informing about the level of information available
#'                       in `votes.election1` and `votes.election2` regarding new entries
#'                       and exits of the election censuses between the two elections or
#'                       indicating how their sum discrepancies should be handled.
#'                       This argument allows the eight options discussed in Pavia (2022)
#'                       as well as an adjusting option. This argument admits nine values: `adjust`,
#'                       `raw`, `regular`, `ordinary`, `simultaneous`, `enriched`, `semifull`,
#'                       `full` and `gold`. See **Details**. Default, `adjust`.
#'
#' @param reference A vector of two components indicating (parties) options in election 1
#'                  and 2, respectively, to be used as reference with `method = "Thomsen"`. This has not effect
#'                  with `method = "IPF"`. The references can be indicated by name or by position.
#'                  If `reference = NULL`, the final solution is constructed as a weighted average
#'                  of all the congruent solutions attained after considering as references all
#'                  combinations of options. Default `NULL`.
#'
#' @param confidence A number between 0 and 1 to be used as level of confidence for the
#'                   confidence intervals of the transition rates. By default `NULL`.
#'                   If `confidence = NULL`, confidence intervals are not computed.
#'
#' @param B An integer indicating the number of samples to be drawn from each crude estimated
#'          confidence interval for estimating final confidence intervals when either
#'          R (J) or C (K) is higher than two. This is not relevant for the 2x2 case.
#'          It can take a while to compute confidence intervals, mainly when `method = Thomsen`.
#'          In general computation burden grows with `B`. Default, `500`.
#'
#' @param Yule.aprox `TRUE`/`FALSE` argument indicating if either Thomsen (1987)'s formula (3.44),
#'                   based on a binormal, or Thomsen (1987)'s formula (3.46), based on Yule's
#'                   approximation of tetrachoric correlation, should be use to estimate
#'                   cross-proportions. Default `FALSE`, formula (3.44).
#'
#' @param tol A number indicating the level of precision to be used to stop the
#'            adjustment of initial/crude count estimates reached using a 2x2 approach in a
#'            general RxC case. This is not relevant for the 2x2 case. Default, `0.000001`.
#'
#' @param ... Other arguments to be passed to the function. Not currently used.
#'
#' @note This function somewhere builds on the .ado (STATA) functions written by Won-ho Park, in 2002.
#'
#' @details Description of the `census.changes` argument in more detail.
#' \itemize{
#'  \item{`adjust`: }{The default value. This is the simplest solution for handling discrepancies
#'                    between the total number of counts for the first and second elections.
#'                    With this value the J1 column-aggregations of the counts
#'                    in `votes.election1` of the first election are proportionally adjusted to
#'                    equal the aggregation of the counts in `votes.election2` of the second election. 
#'                    In this scenario, J is equal to J1 and K equal to K2.}
#'  \item{`raw`: }{This argument accounts for a scenario with two elections elapsed at least
#'                 some months where only the raw election data recorded in the I (territorial) units,
#'                 in which the electoral space under study is divided, are available and net
#'                 entries and net exits are approached from the available information.
#'                 In this scenario, net exits and net entries are estimated according to
#'                 Pavia (2022). When both net entries and exits are no
#'                 null, constraint (15) of Pavia (2022) applies: no transfer between entries and
#'                 exits are allowed. In this scenario, J could be equal to J1 or J1 + 1 and K equal to
#'                 K2 or K2 + 1.}
#'  \item{`simultaneous`: }{This is the value to be used in classical ecological inference problems,
#'                such as in ecological studies of social or racial voting, and in scenarios with two simultaneous elections.
#'                In this scenario, the sum by rows of `votes.election1` and `votes.election2` must coincide.}
#'  \item{`regular`: }{This value accounts for a scenario with
#'                 two elections elapsed at least some months where (i) the column J1
#'                 of `votes.election1` corresponds to new (young) electors who have the right
#'                 to vote for the first time, (ii) net exits and maybe other additional
#'                 net entries are computed according to Pavia (2022). When both net entries and exits
#'                 are no null, constraints (13) and (15) of Pavia (2022) apply. In this scenario, J
#'                 could be equal to J1 or J1 + 1 and K equal to K2 or K2 + 1.}
#'  \item{`ordinary`: }{This value accounts for a scenario
#'                 with two elections elapsed at least some months where (i) the column K1
#'                 of `votes.election2` corresponds to electors who died in the interperiod
#'                 election, (ii) net entries and maybe other additional net exits are
#'                 computed according to Pavia (2022). When both net entries and net exits are no null,
#'                 constraints (14) and (15) of Pavia (2022) apply.
#'                 In this scenario, J could be equal to J1 or J1 + 1 and K equal to K2 or K2 + 1.}
#'  \item{`enriched`: }{This value accounts for a scenario that somewhat combine `regular` and
#'                 `ordinary` scenarios. We consider two elections elapsed at least some months where
#'                 (i) the column J1 of `votes.election1` corresponds to new (young) electors
#'                  who have the right to vote for the first time, (ii) the column K2 of
#'                 `votes.election2` corresponds to electors who died in the interperiod
#'                 election, (iii) other (net) entries and (net) exits are computed according
#'                 to Pavia (2022). When both net entries and net exits are no null, constraints (12) to
#'                 (15) of Pavia (2022) apply. In this scenario, J could be equal
#'                 to J1 or J1 + 1 and K equal to K2 or K2 + 1.}
#'  \item{`semifull`: }{This value accounts for a scenario with two elections elapsed at least some
#'                months, where: (i) the column J1 = J of `votes.election1` totals new
#'                electors (young and immigrants) that have the right to vote for the first time in each polling unit and
#'                (ii) the column K2 = K of `votes.election2` corresponds to total exits of the census
#'                lists (due to death or emigration). In this scenario, the sum by rows of
#'                `votes.election1` and `votes.election2` must agree and constraint (15)
#'                of Pavia (2022) apply.}
#'  \item{`full`: }{This value accounts for a scenario with two elections elapsed at least some
#'                months, where J = J1, K = K2 and (i) the column J - 1 of `votes.election1` totals new (young)
#'                electors that have the right to vote for the first time, (ii) the column J
#'                of `votes.election1` measures new immigrants that have the right to vote and
#'                (iii) the column K of `votes.election2` corresponds to total exits of the census
#'                lists (due to death or emigration). In this scenario, the sum by rows of
#'                `votes.election1` and `votes.election2` must agree and constraints (13)
#'                and (15) of Pavia (2022) apply.}
#'  \item{`gold`: }{This value accounts for a scenario similar to `full`, where J = J1, K = K2 and
#'                 total exits are separated out between exits due to emigration
#'                 (column K - 1 of `votes.election2`) and death (column K of `votes.election2`).
#'                 In this scenario, the sum by rows of `votes.election1` and `votes.election2` must agree.
#'                 Constraints (12) to (15) of Pavia (2022) apply.}
#' }
#'
#' @return
#' A list with the following components
#'  \item{VTM}{ A matrix of order JxK (RxC) with the estimated proportions of the row-standardized vote transitions from election 1 to election 2.
#'              In `raw`, `regular`, `ordinary` and `enriched` scenarios, this matrix includes the row and the column corresponding to net entries
#'              and net exits (when they are present). When `local = TRUE` (default), this matrix is obtained as
#'              composition of the local solutions.}
#'  \item{VTM.votes}{ A matrix of order JxK (RxC) with the estimated vote transfers from election 1 to election 2.
#'                    In `raw`, `regular`, `ordinary` and `enriched` scenarios, this matrix includes the row and the column corresponding to net entries
#'                    and net exits (when they are present). When `local = TRUE` (default), this matrix is obtained as
#'                    aggregation of the local solutions.}
#'  \item{VTM.global}{ A matrix of order JxK (RxC) with the estimated proportions of the row-standardized vote transitions from election 1 to election 2,
#'                     attained directly from the global (whole electoral space) proportions. When `local = FALSE`. `VTM` and `VTM.global` coincide.
#'                     In `raw`, `regular`, `ordinary` and `enriched` scenarios, this matrix includes the row and the column corresponding to net entries
#'                     and net exits (when they are present).}
#'  \item{VTM.votes.global}{ A matrix of order JxK (RxC) with the estimated vote transfers from election 1 to election 2,
#'                           attained directly from the global proportions. When `local = FALSE`, `VTM.votes` and `VTM.votes.global` coincide.
#'                           In `raw`, `regular`, `ordinary` and `enriched` scenarios, this matrix includes the row and the column corresponding to net entries
#'                           and net exits (when they are present).}
#'  \item{VTM.lower}{ A matrix of order JxK (RxC) with the estimated lower limits of the confidence intervals for
#'                      the proportions of the row-standardized vote transitions from election 1 to election 2.
#'                      In `raw`, `regular`, `ordinary` and `enriched` scenarios, this matrix includes the row and the column corresponding to net entries
#'                      and net exits (when they are present). When `confidence = NULL` this is a `NULL` object.}
#'  \item{VTM.upper}{ A matrix of order JxK (RxC) with the estimated upper limits of the confidence intervals for
#'                      the proportions of the row-standardized vote transitions from election 1 to election 2.
#'                      In `raw`, `regular`, `ordinary` and `enriched` scenarios, this matrix includes the row and the column corresponding to net entries
#'                      and net exits (when they are present). When `confidence = NULL` this is a `NULL` object.}
#'  \item{VTM.crude.global}{ A matrix of order JxK (RxC) with the  (inconsistent) crude estimated proportions for the row-standardized
#'                      vote transitions from election 1 to election 2 in the whole space attained in a 2x2 fashion before making them
#'                      consistent using the iterative proportional fitting algorithm or the Thomsen iteratuve algortihm.
#'                      In `raw`, `regular`, `ordinary` and `enriched` scenarios, this matrix includes the row and the 
#'                      column corresponding to net entries and net exits (when they are present).}
#'  \item{VTM.units}{ An array of order JxKxI (RxCxI) with the estimated proportions of the row-standardized vote transitions from election 1 to election 2
#'                    attained for each unit. When `local = FALSE`, this is a `NULL` object.
#'                    In `raw`, `regular`, `ordinary` and `enriched` scenarios, each unit matrix includes the row and the column corresponding to net entries
#'                    and net exits (when they are present).}
#'  \item{VTM.votes.units}{ An array of order JxKxI (RxCxI) with the estimated transfer of votes from election 1 to election 2
#'                    attained for each unit. When `local = FALSE`, this is a `NULL` object.
#'                    In `raw`, `regular`, `ordinary` and `enriched` scenarios, each unit matrix includes the row and the column corresponding to net entries
#'                    and net exits (when they are present).}
#'  \item{VTM.lower.units}{ An array of order JxKxI (RxCxI) with the estimated lower limits of the confidence intervals for
#'                      the proportions of the row-standardized vote transitions from election 1 to election 2 corresponding to each unit.
#'                      When either `local = FALSE` or `confidence = NULL`, this is a `NULL` object.
#'                      In `raw`, `regular`, `ordinary` and `enriched` scenarios, each unit matrix includes the row and the column corresponding to net entries
#'                      and net exits (when they are present).}
#'  \item{VTM.upper.units}{ An array of order JxKxI (RxCxI) with the estimated upper limits of the confidence intervals for
#'                      the proportions of the row-standardized vote transitions from election 1 to election 2 corresponding to each unit.
#'                      When either `local = FALSE` or `confidence = NULL`, this is a `NULL` object.
#'                      In `raw`, `regular`, `ordinary` and `enriched` scenarios, each unit matrix includes the row and the column corresponding to net entries
#'                      and net exits (when they are present).}
#'  \item{VTM.crude.units}{ An array of order JxKxI (RxCxI) with the (inconsistent) crude estimated proportions of the row-standardized vote transitions from election 1 to election 2
#'                    attained for each unit in a 2x2 fashion before making them consistent using the iterative proportional fitting algorithm or the Thomsen iterative algorithm. 
#'                    When `local = FALSE`, this is a NULL object. In `raw`, `regular`, `ordinary` and `enriched` scenarios, each unit matrix includes the row and the column 
#'                    corresponding to net entries and net exits (when they are present).}
#'  \item{correlations}{ A matrix of order JxK (Rxc) with the across units correlations between options for the proportions
#'                       in the transformed scale.}
#'  \item{reference.outputs}{ A list with three components: `vjk.averages`, `vjk.by.reference` and `vjk.units.by.reference`.
#'                            The first component `vjk.averages` is a JxKx8 array with eight different global solutions
#'                            of transfer matrix of votes attained after combining with different weights each of the
#'                            solutions obtained using the different combinations of a row and a column option as reference.
#'                            The second component `vjk.by.reference` is a JxKx(J1·K1) array with the J1K1 different
#'                            global solutions of transfer matrix of votes attained after choosing as reference all the
#'                            possible combination of a row and a column option. The third component `vjk.units.by.reference`
#'                            is a JxKxIx(J1·K1) array with the local solutions linked to `vjk.by.reference`.
#'                            When either `method = "IPF"` or `reference` is not `NULL`, this is a `NULL` object.}
#'
#'  \item{iter}{ A vector of either length 1 (when `reference` is different of `NULL`) or J1·K1 with
#'              the number of iterations needed by the Thomsen algorithm to reach convergence for
#'              each reference pair. When `method = "Thomsen"` this is a `NULL` object.}
#'
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'
#' @export
#'
#' @family latent structure ecological inference functions
#'
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
#' example <- ecolRxC(votes1, votes2, method = "IPF")$VTM
#'
#' @importFrom stats pnorm qnorm quantile runif


ecolRxC <- function(votes.election1,
                    votes.election2,
                    scale = "probit",
                    method = "Thomsen",
                    local = TRUE,
                    census.changes = c("adjust", "raw", "regular", "ordinary", "enriched",
                                       "simultaneous", "semifull", "full", "gold"),
                    reference = NULL,
                    confidence = NULL,
                    B = 500,
                    Yule.aprox = FALSE,
                    tol = 0.000001,
                    ...){

  argg <- c(as.list(environment()), list(...))
  nothing <- tests_inputs_ecolRxC(argg)

  if(method == "IPF"){
    if(local){
      output <- ecolRxC_local(votes.election1 = votes.election1,
                              votes.election2 = votes.election2,
                              scale = scale,
                              census.changes = census.changes,
                              confidence = confidence,
                              B = B,
                              Yule.aprox = Yule.aprox,
                              tol = tol,
                              ...)
    } else {
      output <- ecolRxC_basic(votes.election1 = votes.election1,
                              votes.election2 = votes.election2,
                              scale = scale,
                              census.changes = census.changes,
                              confidence = confidence,
                              B = B,
                              Yule.aprox = Yule.aprox,
                              tol = tol,
                              ...)

    }
  } else {
    output <- ecolRxC_Thomsen(votes.election1 = votes.election1,
                              votes.election2 = votes.election2,
                              scale = scale,
                              census.changes = census.changes,
                              reference = reference,
                              confidence = confidence,
                              B = B,
                              Yule.aprox = Yule.aprox,
                              tol = tol,
                              ...)

  }

  class(output) <- "ecolRxC"
  return(output)
}





ecolRxC_basic <- function(votes.election1,
                          votes.election2,
                          scale = "logit",
                          census.changes = c("adjust", "raw", "regular", "ordinary", "enriched",
                                          "simultaneous", "semifull", "full", "gold"),
                          confidence = 0.95,
                          B = 500,
                          Yule.aprox = FALSE,
                          tol = 0.000001,
                          ...){

  # argg <- c(as.list(environment()), list(...))
  x0 <- as.matrix(votes.election1)
  y0 <- as.matrix(votes.election2)


  # inputs
  inputs <- list("votes.election1" = votes.election1, "votes.election2" = votes.election2,
                 "scale" = scale, method = "IPF", "local" = FALSE, "census.changes" = census.changes[1],
                 "confidence" = confidence, "B" = B, "Yule.aprox" = Yule.aprox, "tol" = tol)

  census.changes <- scenario <- census.changes[1]

  # Basic 2x2 function
  if (scale == "logit"){
     ecol2x2 <- Thomsen_logit_2x2
  } else {
     ecol2x2 <- Thomsen_probit_2x2
  }

  # Data preparation
  net <- compute_net_voters(x0 = x0, y0 = y0, scenario = scenario)
  x <- net$x
  y <- net$y
  vector.columna <- colSums(y)
  vector.fila <- colSums(x)

  # Parameters
  J <- ncol(x)
  K <- ncol(y)
  J0 <- ncol(x0)
  K0 <- ncol(y0)

  # Names of election options
  names1 <- colnames(x)
  names2 <- colnames(y)

  # Estimation of crude transitions
  if (J == 2L & K == 2L){
    mt <- ecol2x2(x[, 1L], rowSums(x), y[, 1L], rowSums(y),
                  confidence = confidence, Yule.aprox = Yule.aprox)
    correlations <- matrix(c(mt$cor, -mt$cor), 2L, 2L)
    VTM.crude <- mt$PTM
    VTM.crude.l <- mt$PTM_low
    VTM.crude.u <- mt$PTM_high
    VTM.votes <- VTM.crude * colSums(x)
    if (is.null(confidence)){
      VTM.crude.l <- VTM.crude.u <- NULL
    }

    output <- list("VTM" = VTM.crude, "VTM.votes" = VTM.votes, "VTM.global" = VTM.crude,
                   "VTM.votes.global" = VTM.crude,
                   "VTM.lower" = VTM.crude.l, "VTM.upper" = VTM.crude.u, "VTM.crude.global" = VTM.crude,
                   "VTM.units" = NULL, "VTM.votes.units" = NULL, "VTM.lower.units" = NULL,
                   "VTM.upper.units" = NULL, "VTM.crude.units" = NULL, "correlations" = correlations,
                   "reference.outputs" = NULL, "iter" = NULL, "inputs" = inputs)
    return(output)
  } else {
    correlations <- VTM.crude <- VTM.crude.l <- VTM.crude.u <- matrix(NA, J, K)
    for (j in 1L:J){
      for (k in 1L:K){
        mt <- ecol2x2(x[, j], rowSums(x), y[, k], rowSums(y),
                      confidence = confidence, Yule.aprox = Yule.aprox)
        VTM.crude[j, k] <- mt$PTM[1L, 1L]
        VTM.crude.l[j, k] <- mt$PTM_low[1L, 1L]
        VTM.crude.u[j, k] <- mt$PTM_high[1L, 1L]
        correlations[j, k] <- mt$cor
      }
    }
  }

  # Constraints given by the scenario
  if (census.changes == "adjust") {
    vector.fila <- vector.fila * sum(vector.columna)/sum(vector.fila)
  } else {
    mt <- constraints_scenario(VTM.crude, VTM.crude.l, VTM.crude.u, scenario, J0, K0)
    VTM.crude <- mt$VTM.crude
    VTM.crude.l <- mt$VTM.crude.l
    VTM.crude.u <- mt$VTM.crude.u
  }
  # Adjustment of initial solutions for congruence
  VTM.votes <- IPF(VTM.crude, vector.columna, vector.fila, precision = tol)
  VTM <- VTM.votes/rowSums(VTM.votes)

  # Uncertainty
  VTM.lower <- VTM.upper <- NULL
  if(!is.null(confidence)){
    muestra <- extract_sample(TM.low = VTM.crude.l, TM.upp = VTM.crude.u, B = B)
    intervalos <- interval_transitions(muestra = muestra,
                                       vector.fila = vector.fila,
                                       vector.columna = vector.columna,
                                       tol = tol,
                                       confidence = confidence)
    VTM.lower <- intervalos$TM.low
    VTM.upper <- intervalos$TM.upp
    colnames(VTM.lower) <- colnames(VTM.upper) <- names2
    rownames(VTM.lower) <- rownames(VTM.upper) <- names1
  }

  # Outputs
  colnames(VTM) <- colnames(VTM.votes) <- colnames(VTM.crude) <- colnames(correlations) <- names2
  rownames(VTM) <- rownames(VTM.votes) <- rownames(VTM.crude) <- rownames(correlations) <- names1

  output <- list("VTM" = VTM, "VTM.votes" = VTM.votes, "VTM.global" = VTM, "VTM.votes.global" = VTM.votes,
                 "VTM.lower" = VTM.lower, "VTM.upper" = VTM.upper, "VTM.crude.global" = VTM.crude,
                 "VTM.units" = NULL, "VTM.votes.units" = NULL, "VTM.lower.units" = NULL,
                 "VTM.upper.units" = NULL, "VTM.crude.units" = NULL, "correlations" = correlations,
                 "reference.outputs" = NULL, "iter" = NULL, "inputs" = inputs)
  return(output)
}


ecolRxC_local <- function(votes.election1,
                          votes.election2,
                          scale = "logit",
                          census.changes = c("adjust", "raw", "regular", "ordinary", "enriched",
                                             "simultaneous", "semifull", "full", "gold"),
                          confidence = 0.95,
                          B = 500,
                          Yule.aprox = FALSE,
                          tol = 0.000001,
                          ...){

  # argg <- c(as.list(environment()), list(...))
  x0 <- as.matrix(votes.election1)
  y0 <- as.matrix(votes.election2)

  # inputs
  inputs <- list("votes.election1" = votes.election1, "votes.election2" = votes.election2,
                 "scale" = scale, method = "IPF", "local" = TRUE, "census.changes" = census.changes[1],
                 "confidence" = confidence, "B" = B, "Yule.aprox" = Yule.aprox, "tol" = tol)

  census.changes <- scenario <- census.changes[1]

  # Basic 2x2 function
  if (scale == "logit"){
    ecol2x2 <- Thomsen_local_logit_2x2
  } else {
    ecol2x2 <- Thomsen_local_probit_2x2
  }

  # Data preparation
  net <- compute_net_voters(x0 = x0, y0 = y0, scenario = scenario)
  x <- net$x
  y <- net$y
  vector.columna <- colSums(y)
  vector.fila <- colSums(x)

  # Parameters
  I <- nrow(x)
  J <- ncol(x)
  K <- ncol(y)
  J0 <- ncol(x0)
  K0 <- ncol(y0)

  # Names of election options and units
  names1 <- colnames(x)
  names2 <- colnames(y)
  names.units <- rownames(x)
  nombres <- c(list(names1), list(names2), list(names.units))

  # Estimation of crude transitions
  if (J == 2L & K == 2L){
    mt <- ecol2x2(x[, 1L], rowSums(x), y[, 1L], rowSums(y),
                  confidence = confidence, Yule.aprox = Yule.aprox)
    correlations <- matrix(c(mt$cor, -mt$cor), 2L, 2L)
    VTM.crude <- mt$PTM
    VTM.votes.global <- VTM.crude * colSums(x)
    VTM.crude.l <- mt$PTM_low
    VTM.crude.u <- mt$PTM_high
    VTM.crude.local <- mt$PTM_local
    VTM.crude.l.local <- mt$PTM_low_local
    VTM.crude.u.local <- mt$PTM_high_local
    VTM.votes.units <- VTM.crude.local * array(as.vector(kronecker(rep(1L, 2L), t(x))),
                                               dim(VTM.crude.local))
    VTM.votes <- apply(VTM.votes.units, c(1L, 2L), sum)
    VTM <- VTM.votes/colSums(x)
    VTM.votes.l.local <- VTM.crude.l.local * array(as.vector(kronecker(rep(1L, 2L), t(x))),
                                                   dim(VTM.crude.local))
    VTM.lower <- apply(VTM.votes.l.local, c(1L, 2L), sum)/colSums(x)
    VTM.votes.u.local <- VTM.crude.u.local * array(as.vector(kronecker(rep(1L, 2L), t(x))),
                                                   dim(VTM.crude.local))
    VTM.upper <- apply(VTM.votes.u.local, c(1L, 2L), sum)/colSums(x)
    if (is.null(confidence)){
      VTM.lower <- VTM.upper <- VTM.crude.l.local <- VTM.crude.u.local <- NULL
    }
    output <- list("VTM" = VTM, "VTM.votes" = VTM.votes, "VTM.global" = VTM.crude,
                   "VTM.votes.global" = VTM.votes.global,
                   "VTM.lower" = VTM.lower, "VTM.upper" = VTM.upper, "VTM.crude.global" = VTM.crude,
                   "VTM.units" = VTM.crude.local, "VTM.votes.units" = VTM.votes.units,
                   "VTM.lower.units" = VTM.crude.l.local, "VTM.upper.units" = VTM.crude.u.local,
                   "VTM.crude.units" = VTM.crude.local, "correlations" = correlations,
                   "reference.outputs" = NULL, "iter" = NULL,
                   "inputs" = inputs)
    return(output)
  } else {
    correlations <- VTM.crude <- VTM.crude.l <- VTM.crude.u <- matrix(NA, J, K)
    VTM.crude.local <- VTM.crude.l.local <- VTM.crude.u.local <- array(NA, c(J, K, I))
    for (j in 1L:J){
      for (k in 1L:K){
        mt <- ecol2x2(x[, j], rowSums(x), y[, k], rowSums(y),
                      confidence = confidence, Yule.aprox = Yule.aprox)
        VTM.crude[j, k] <- mt$PTM[1L, 1L]
        VTM.crude.l[j, k] <- mt$PTM_low[1L, 1L]
        VTM.crude.u[j, k] <- mt$PTM_high[1L, 1L]
        correlations[j, k] <- mt$cor
        VTM.crude.local[j, k, ] <- mt$PTM_local[1L, 1L, ]
        VTM.crude.l.local[j, k, ] <- mt$PTM_low_local[1L, 1L, ]
        VTM.crude.u.local[j, k, ] <- mt$PTM_high_local[1L, 1L, ]
      }
    }
  }

  # Constraints given by the scenario
  if (census.changes == "adjust") {
    vector.fila <- vector.fila * sum(vector.columna)/sum(vector.fila)
    x <- t(t(x) * rowSums(y)/rowSums(x))
  } else {
    mt <- constraints_scenario(VTM.crude, VTM.crude.l, VTM.crude.u, scenario, J0, K0)
    VTM.crude <- mt$VTM.crude
    VTM.crude.l <- mt$VTM.crude.l
    VTM.crude.u <- mt$VTM.crude.u
    mt <- constraints_scenario_local(VTM.crude.local, VTM.crude.l.local,
                                     VTM.crude.u.local, scenario, J0, K0)
    VTM.crude.local <- mt$VTM.crude.local
    VTM.crude.l.local <- mt$VTM.crude.l.local
    VTM.crude.u.local <- mt$VTM.crude.u.local
  }

  # Adjustment of initial solutions for congruence
  # Global matrix
  VTM.votes.global <- IPF(VTM.crude, vector.columna, vector.fila, precision = tol)
  VTM.global <- VTM.votes.global/rowSums(VTM.votes.global)
  # Units
  VTM.votes.units <- VTM.units <- VTM.crude.local
  for (ii in 1L:I){
    VTM.votes.units[, , ii] <- IPF(VTM.crude.local[, , ii], y[ii, ], x[ii, ], precision = tol)
    VTM.units[, , ii] <- VTM.votes.units[, , ii]/rowSums(VTM.votes.units[, , ii])
    VTM.units[x[ii, ] == 0L, , ii] <- 0L
  }
  # Composition/aggregation matrix
  VTM.votes <- apply(VTM.votes.units, c(1L, 2L), sum)
  VTM <- VTM.votes/rowSums(VTM.votes)

  # Uncertainty
  VTM.lower <- VTM.upper <- NULL
  if(!is.null(confidence)){
    VTM.lower <- VTM.upper <- matrix(0L, J, K)
    for (ii in 1L:I){
      muestra <- extract_sample(TM.low = VTM.crude.l.local[, , ii],
                                TM.upp = VTM.crude.u.local[, , ii],
                                B = B)
      intervalos <- interval_transfers(muestra = muestra,
                                      vector.fila = x[ii, ],
                                      vector.columna = y[ii, ],
                                      tol = tol,
                                      confidence = confidence)
      VTM.crude.l.local[, , ii] <- intervalos$TM.low/rowSums(VTM.votes.units[, , ii])
      VTM.crude.u.local[, , ii] <- intervalos$TM.upp/rowSums(VTM.votes.units[, , ii])
      VTM.crude.l.local[x[ii, ] == 0L, , ii] <- 0L
      VTM.crude.u.local[x[ii, ] == 0L, , ii] <- 0L
      VTM.lower <- VTM.lower + intervalos$TM.low
      VTM.upper <- VTM.upper + intervalos$TM.upp
    }
    VTM.lower <- VTM.lower/rowSums(VTM.votes)
    VTM.upper <- VTM.upper/rowSums(VTM.votes)
    colnames(VTM.lower) <- colnames(VTM.upper) <- names2
    rownames(VTM.lower) <- rownames(VTM.upper) <- names1
    dimnames(VTM.crude.l.local) <- dimnames(VTM.crude.u.local) <- nombres
  }

  # Outputs
  colnames(VTM) <- colnames(VTM.votes) <- colnames(VTM.global) <- colnames(VTM.votes.global) <-
    colnames(VTM.crude) <- colnames(correlations) <- names2
  rownames(VTM) <- rownames(VTM.votes) <- rownames(VTM.global) <- rownames(VTM.votes.global) <-
     rownames(VTM.crude) <- rownames(correlations) <- names1

  dimnames(VTM.units) <- dimnames(VTM.votes.units) <- dimnames(VTM.crude.local) <- nombres

  if (is.null(confidence)){
    VTM.lower <- VTM.upper <- VTM.crude.l.local <- VTM.crude.u.local <- NULL
  }

  output <- list("VTM" = VTM, "VTM.votes" = VTM.votes, "VTM.global" = VTM.global,
                 "VTM.votes.global" = VTM.votes.global,
                 "VTM.lower" = VTM.lower, "VTM.upper" = VTM.upper, "VTM.crude.global" = VTM.crude,
                 "VTM.units" = VTM.units, "VTM.votes.units" = VTM.votes.units,
                 "VTM.lower.units" = VTM.crude.l.local, "VTM.upper.units" = VTM.crude.u.local,
                 "VTM.crude.units" = VTM.crude.local, "correlations" = correlations,
                 "reference.outputs" = NULL, "iter" = NULL, "inputs" = inputs)

  return(output)
}



ecolRxC_Thomsen <- function(votes.election1,
                            votes.election2,
                            scale = "logit",
                            census.changes = c("adjust", "raw", "regular", "ordinary", "enriched",
                                               "simultaneous", "semifull", "full", "gold"),
                            reference = NULL,
                            confidence = 0.95,
                            B = 500,
                            Yule.aprox = FALSE,
                            tol = 0.000001,
                            ...){

  # reference: A vector of two components indicating the reference (parties) options in election 1 and 2, respectively.
  #            The references can be indicated by name or by position.
  #            If NULL, the solution is constructed as a mean of all the solutions attained after considering
  #            as references all combinations of options.

  # argg <- c(as.list(environment()), list(...))
  x0 <- as.matrix(votes.election1)
  y0 <- as.matrix(votes.election2)

  # inputs
  inputs <- list("votes.election1" = votes.election1, "votes.election2" = votes.election2,
                 "reference" = reference, "scale" = scale, method = "Thomsen", "local" = TRUE,
                 "census.changes" = census.changes[1], "confidence" = confidence,
                 "B" = B, "Yule.aprox" = Yule.aprox, "tol" = tol)

  census.changes <- scenario <- census.changes[1]

  # Basic 2x2 function
  if (scale == "logit"){
    ecol2x2 <- Thomsen_local_logit_2x2
  } else {
    ecol2x2 <- Thomsen_local_probit_2x2
  }

  # Data preparation
  net <- compute_net_voters(x0 = x0, y0 = y0, scenario = scenario)
  x <- net$x
  y <- net$y
  vector.columna <- colSums(y)
  vector.fila <- colSums(x)

  # Parameters
  I <- nrow(x)
  J <- ncol(x)
  K <- ncol(y)
  J0 <- ncol(x0)
  K0 <- ncol(y0)

  # Names of election options and units
  names1 <- colnames(x)
  names2 <- colnames(y)
  names.units <- rownames(x)

  # Estimation of crude cross-probabilities
  if (J == 2L & K == 2L){
    mt <- ecol2x2(x[, 1L], rowSums(x), y[, 1L], rowSums(y),
                  confidence = confidence, Yule.aprox = Yule.aprox)
    correlations <- matrix(c(mt$cor, -mt$cor), 2L, 2L)
    VTM.crude <- mt$PTM
    VTM.votes.global <- VTM.crude * colSums(x)
    VTM.crude.l <- mt$PTM_low
    VTM.crude.u <- mt$PTM_high
    VTM.crude.local <- mt$PTM_local
    VTM.crude.l.local <- mt$PTM_low_local
    VTM.crude.u.local <- mt$PTM_high_local
    VTM.votes.units <- VTM.crude.local * array(as.vector(kronecker(rep(1L, 2L), t(x))),
                                               dim(VTM.crude.local))
    VTM.votes <- apply(VTM.votes.units, c(1L, 2L), sum)
    VTM <- VTM.votes/colSums(x)
    VTM.votes.l.local <- VTM.crude.l.local * array(as.vector(kronecker(rep(1L, 2L), t(x))),
                                                   dim(VTM.crude.local))
    VTM.lower <- apply(VTM.votes.l.local, c(1L, 2L), sum)/colSums(x)
    VTM.votes.u.local <- VTM.crude.u.local * array(as.vector(kronecker(rep(1L, 2L), t(x))),
                                                   dim(VTM.crude.local))
    VTM.upper <- apply(VTM.votes.u.local, c(1L, 2L), sum)/colSums(x)
    if (is.null(confidence)){
      VTM.lower <- VTM.upper <- VTM.crude.l.local <- VTM.crude.u.local <- NULL
    }
    output <- list("VTM" = VTM, "VTM.votes" = VTM.votes, "VTM.global" = VTM.crude,
                   "VTM.votes.global" = VTM.votes.global,
                   "VTM.lower" = VTM.lower, "VTM.upper" = VTM.upper, "VTM.crude.global" = VTM.crude,
                   "VTM.units" = VTM.crude.local, "VTM.votes.units" = VTM.votes.units,
                   "VTM.lower.units" = VTM.crude.l.local, "VTM.upper.units" = VTM.crude.u.local,
                   "VTM.crude.units" = VTM.crude.local, "correlations" = correlations,
                   "reference.outputs" = NULL, "iter" = NULL,
                   "inputs" = inputs)
    return(output)
  } else {
    correlations <- VTM.crude <- VTM.crude.l <- VTM.crude.u <- matrix(NA, J, K)
    VTM.crude.local <- pjk.crude.local <- pjk.crude.l.local <- pjk.crude.u.local <- array(NA, c(J, K, I))
    for (j in 1L:J){
      for (k in 1L:K){
        mt <- ecol2x2(x[, j], rowSums(x), y[, k], rowSums(y),
                      confidence = confidence, Yule.aprox = Yule.aprox)
        VTM.crude[j, k] <- mt$PTM[1L, 1L]
        VTM.crude.l[j, k] <- mt$PTM_low[1L, 1L]
        VTM.crude.u[j, k] <- mt$PTM_high[1L, 1L]
        correlations[j, k] <- mt$cor
        VTM.crude.local[j, k, ] <- mt$PTM_local[1L, 1L, ]
        pjk.crude.local[j, k, ] <- mt$pjk[1L, 1L, ]
        pjk.crude.l.local[j, k, ] <- mt$pjk_low[1L, 1L, ]
        pjk.crude.u.local[j, k, ] <- mt$pjk_high[1L, 1L, ]
      }
    }
  }

  # Constraints given by the scenario
  if (census.changes != "adjust") {
    mt <- constraints_scenario(VTM.crude, VTM.crude.l, VTM.crude.u, scenario, J0, K0)
    VTM.crude <- mt$VTM.crude
  }
  colnames(VTM.crude) <- names2
  rownames(VTM.crude) <- names1

  # Null constraints given by the scenario and by aggregate null rows or columns
  mt <- constraints_zeros_local(pjk.crude.local, pjk.crude.l.local, pjk.crude.u.local,
                                scenario, J0, K0, x, y)
  pjk.crude.local <- mt$pjk.crude.local
  pjk.crude.l.local <- mt$pjk.crude.l.local
  pjk.crude.u.local <- mt$pjk.crude.u.local

  # Adjustment of initial solutions for congruence according to the Thomsen algorithm
  T.adj <- Thomsen_iter_algorithm(pjk.crude.local = pjk.crude.local, Yule.aprox = Yule.aprox,
                                  reference = reference, scale = scale, x = x, y = y,
                                  J0 = J0, K0 = K0, tol = tol)
  vjk.units.multi <- T.adj$vjk.units.multi
  iter <- T.adj$iter
  vjk.units <- T.adj$vjk.units

  if (is.null(reference)){
    # Solutions combining reference solutions using abs(correlations) as weights
    weights <- abs(as.vector(t(correlations[1L:J0, 1L:K0])))
    weights <- weights/sum(weights)
    W <- array(rep(weights, each = I*J*K), dim(vjk.units.multi))
    vjk.units <- apply(vjk.units.multi * W, c(1L, 2L, 3L), sum)
  }

  # Units
  VTM.units <- vjk.units
  for (ii in 1L:I){
    VTM.units[, , ii] <- vjk.units[, , ii]/rowSums(vjk.units[, , ii])
    VTM.units[x[ii, ] == 0L, , ii] <- 0L
  }

  # Global matrix
  # vjk_r1_r2 <- solutions_by_reference(vjk.units.multi, x0, y0)
  # W <- array(rep(weights, each = J0*K0), dim(vjk_r1_r2))
  # VTM.votes <- VTM.votes.global <- apply(vjk_r1_r2 * W, c(1L, 2L), sum)
  VTM.votes <- VTM.votes.global <- apply(vjk.units, c(1L, 2L), sum)
  VTM <- VTM.global <- VTM.votes.global/rowSums(VTM.votes.global)

  # Reference outputs (other solutions)
  reference.outputs <- NULL
  if (is.null(reference)){
    reference.outputs <- reference_outputs(vjk.units.multi = vjk.units.multi, VTM.crude = VTM.crude,
                                           pjk.crude.local = pjk.crude.local, x0 = x0, y0 = y0,
                                           correlations = correlations)
  }

  # Uncertainty
  VTM.upper <- VTM.lower <- matrix(NA, J, K)
  if(!is.null(confidence)){
    muestra_vjk_local <- array(NA, c(J, K, I, B))
    for (bb in 1L:B){
      muestra_vjk_local[, , , bb] <- extract_muestra_vjk(pjk.low = pjk.crude.l.local,
                                                         pjk.upp = pjk.crude.u.local,
                                                         Yule.aprox = Yule.aprox,
                                                         reference = reference,
                                                         scale = scale, x = x, y = y,
                                                         J0 = J0, K0 = K0, tol = tol,
                                                         correlations = correlations)
    }

    # Confidence Intervals
    int.conf <- intervals_Thomsen(muestra_vjk_local = muestra_vjk_local, vjk.units = vjk.units,
                                  x = x, confidence = confidence)
    pjk.crude.l.local <- int.conf$VTM.l.local
    pjk.crude.u.local <- int.conf$VTM.u.local
    VTM.lower <- int.conf$VTM.lower
    VTM.upper <- int.conf$VTM.upper
  }

  # Outputs
  colnames(VTM) <- colnames(VTM.votes) <- colnames(VTM.lower) <-
    colnames(VTM.global) <- colnames(VTM.votes.global) <-
    colnames(VTM.upper) <- colnames(VTM.crude) <- colnames(correlations) <- names2
  rownames(VTM) <- rownames(VTM.votes) <- rownames(VTM.lower) <-
    rownames(VTM.global) <- rownames(VTM.votes.global) <-
    rownames(VTM.upper) <- rownames(VTM.crude) <- rownames(correlations) <- names1

  nombres <- c(dimnames(VTM), list(names.units))
  dimnames(pjk.crude.l.local) <- dimnames(pjk.crude.u.local) <-
    dimnames(VTM.units) <- dimnames(vjk.units) <- dimnames(VTM.crude.local) <- nombres

  if (is.null(confidence)){
    VTM.lower <- VTM.upper <- pjk.crude.l.local <- pjk.crude.u.local <- NULL
  }

  output <- list("VTM" = VTM, "VTM.votes" = VTM.votes, "VTM.global" = VTM.global,
                 "VTM.votes.global" = VTM.votes.global,
                 "VTM.lower" = VTM.lower, "VTM.upper" = VTM.upper, "VTM.crude.global" = VTM.crude,
                 "VTM.units" = VTM.units, "VTM.votes.units" = vjk.units,
                 "VTM.lower.units" = pjk.crude.l.local, "VTM.upper.units" = pjk.crude.u.local,
                 "VTM.crude.units" = VTM.crude.local, "correlations" = correlations,
                 "reference.outputs" = reference.outputs, "iter" = iter, "inputs" = inputs)

  return(output)
}


