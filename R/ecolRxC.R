#' Ecological Inference of RxC Tables by Latent Structure Approaches
#'
#' @description  Estimates JxK (RxC) vote transfer matrices (ecological contingency tables) based on Thomsen (1987) and Park (2008) approaches.
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @references Achen, C.H. (2000). The Thomsen Estimator for Ecological Inference (Unpublished manuscript). University of Michigan.
#' @references Park, W.-H. (2008). Ecological Inference and Aggregate Analysis of Elections. PhD Dissertation. University of Michigan.
#' @references Pavía, J.M. (2023). Adjustment of initial estimates of voter transition probabilities to guarantee consistency and completeness. *SN Social Sciences*, 3, 75. \doi{10.1007/s43545-023-00658-y}.
#' @references Pavía, J.M. and Thomsen, S.R. (2025). ecolRxC: Ecological inference estimation of RxC tables using latent structure approaches. *Political Science Research and Methods*, 13(4), 943-961. \doi{10.1017/psrm.2024.57}.
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
#'                       as well as an adjusting option. This argument admits nine values: `"adjust"`,
#'                       `"raw"`, `"regular"`, `"ordinary"`, `"simultaneous"`, `"enriched"`, `"semifull"`,
#'                       `"full"` and `"gold"`. See **Details**. Default, `"adjust"`.
#'
#' @param reference A vector of two components indicating (parties) options in election 1
#'                  and 2, respectively, to be used as reference with `method = "Thomsen"`. This has not effect
#'                  with `method = "IPF"`. The references can be indicated by name or by position.
#'                  If `reference = NULL`, the final solution is constructed as a weighted average
#'                  of all the congruent solutions attained after considering as references all
#'                  combinations of options. Default `NULL`.
#'                  
#' @param regions Optional argument used to partition the set of units into subsets. It must be a vector 
#'                of length equal to the number of units, indicating the subset (region) to which each 
#'                unit belongs. When provided, the selected ecolRxC specification is then applied independently 
#'                within each subset and the joint estimates are obtained by properly combining/aggregating the 
#'                different region estimates. By default `NULL`, no partition.
#'
#' @param confidence A number between 0 and 1 to be used as level of confidence for the
#'                   confidence intervals of the transition rates. By default `NULL`.
#'                   If `confidence = NULL`, confidence intervals are not computed.
#'
#' @param B Either a non-positive number or a positive integer number. When `B` is non-positive,
#'          confidence intervals are computed in a similar way than in Pavia-Miralles (2005),
#'          <https://www.jstor.org/stable/27590658>, with the methodology adapted to this problem.
#'          When `B` is a positive integer number, it represents the number of samples to 
#'          be drawn from each crude estimated confidence interval for estimating final 
#'          confidence intervals when either R (J) or C (K) is higher than two. 
#'          This argument is not relevant for the 2x2 case.
#'          It can take a while to compute confidence intervals, mainly when `method = Thomsen`
#'          and `B` is a large positive integer; in general computation burden grows as 
#'          a function of `B`, the number of units and the dimension of the transition table. 
#'          Default, `0`.
#'
#' @param Yule.aprox `TRUE`/`FALSE` argument indicating if either Thomsen (1987)'s formula (3.44),
#'                   based on a binormal, or Thomsen (1987)'s formula (3.46), based on Yule's
#'                   approximation of tetrachoric correlation, should be use to estimate
#'                   cross-proportions. Default `FALSE`, formula (3.44).
#'
#' @param ref.combination A character string specifying how unit table estimates should be
#'                        combined to generate the default solution when `method = "Thomsen"` is used
#'                        and the estimates are obtained using all cells as reference.
#'                        This argument allows for the eight options discussed in Pavia and Thomsen (2025):
#'                        `"Mean"`, `"RCNV"`, `"SQRCNV"`, `"SQRM"`, `"AVCR"`, `"LRCNV"`, `"LSQRCNV"`,
#'                        and `"LSQRM"`; \doi{10.1017/psrm.2024.57}. The default is `"LRCNV"`. 
#'                        The eight solutions are included in the component `vjk.averages` of `reference.outputs`.
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
#' \describe{
#'  \item{`adjust`:}{ The default value. This is the simplest solution for handling discrepancies
#'                    between the total number of counts for the first and second elections.
#'                    With this value the J1 column-aggregations of the counts
#'                    in `votes.election1` of the first election are proportionally adjusted to
#'                    equal the aggregation of the counts in `votes.election2` of the second election. 
#'                    In this scenario, J is equal to J1 and K equal to K2.}
#'  \item{`raw`:}{ This argument accounts for a scenario with two elections elapsed at least
#'                 some months where only the raw election data recorded in the I (territorial) units,
#'                 in which the electoral space under study is divided, are available and net
#'                 entries and net exits are approached from the available information.
#'                 In this scenario, net exits and net entries are estimated according to
#'                 Pavia (2022). When both net entries and exits are no
#'                 null, constraint (15) of Pavia (2022) applies: no transfer between entries and
#'                 exits are allowed. In this scenario, J could be equal to J1 or J1 + 1 and K equal to
#'                 K2 or K2 + 1.}
#'  \item{`simultaneous`:}{ This is the value to be used in classical ecological inference problems,
#'                such as in ecological studies of social or racial voting, and in scenarios with two simultaneous elections.
#'                In this scenario, the sum by rows of `votes.election1` and `votes.election2` must coincide.}
#'  \item{`regular`:}{ This value accounts for a scenario with
#'                 two elections elapsed at least some months where (i) the column J1
#'                 of `votes.election1` corresponds to new (young) electors who have the right
#'                 to vote for the first time, (ii) net exits and maybe other additional
#'                 net entries are computed according to Pavia (2022). When both net entries and exits
#'                 are no null, constraints (13) and (15) of Pavia (2022) apply. In this scenario, J
#'                 could be equal to J1 or J1 + 1 and K equal to K2 or K2 + 1.}
#'  \item{`ordinary`:}{ This value accounts for a scenario
#'                 with two elections elapsed at least some months where (i) the column K1
#'                 of `votes.election2` corresponds to electors who died in the interperiod
#'                 election, (ii) net entries and maybe other additional net exits are
#'                 computed according to Pavia (2022). When both net entries and net exits are no null,
#'                 constraints (14) and (15) of Pavia (2022) apply.
#'                 In this scenario, J could be equal to J1 or J1 + 1 and K equal to K2 or K2 + 1.}
#'  \item{`enriched`:}{ This value accounts for a scenario that somewhat combine `regular` and
#'                 `ordinary` scenarios. We consider two elections elapsed at least some months where
#'                 (i) the column J1 of `votes.election1` corresponds to new (young) electors
#'                  who have the right to vote for the first time, (ii) the column K2 of
#'                 `votes.election2` corresponds to electors who died in the interperiod
#'                 election, (iii) other (net) entries and (net) exits are computed according
#'                 to Pavia (2022). When both net entries and net exits are no null, constraints (12) to
#'                 (15) of Pavia (2022) apply. In this scenario, J could be equal
#'                 to J1 or J1 + 1 and K equal to K2 or K2 + 1.}
#'  \item{`semifull`:}{ This value accounts for a scenario with two elections elapsed at least some
#'                months, where: (i) the column J1 = J of `votes.election1` totals new
#'                electors (young and immigrants) that have the right to vote for the first time in each polling unit and
#'                (ii) the column K2 = K of `votes.election2` corresponds to total exits of the census
#'                lists (due to death or emigration). In this scenario, the sum by rows of
#'                `votes.election1` and `votes.election2` must agree and constraint (15)
#'                of Pavia (2022) apply.}
#'  \item{`full`:}{ This value accounts for a scenario with two elections elapsed at least some
#'                months, where J = J1, K = K2 and (i) the column J - 1 of `votes.election1` totals new (young)
#'                electors that have the right to vote for the first time, (ii) the column J
#'                of `votes.election1` measures new immigrants that have the right to vote and
#'                (iii) the column K of `votes.election2` corresponds to total exits of the census
#'                lists (due to death or emigration). In this scenario, the sum by rows of
#'                `votes.election1` and `votes.election2` must agree and constraints (13)
#'                and (15) of Pavia (2022) apply.}
#'  \item{`gold`:}{ This value accounts for a scenario similar to `full`, where J = J1, K = K2 
#'                 where (i) the column J - 1 of `votes_election1` totals new young
#'                 electors that have the right to vote for the first time, (ii) the column J
#'                 of `votes_election1` measures new immigrants that have the right to vote,
#'                 and total exits are separated out between (iii) exits due to emigration
#'                 (column K - 1 of `votes.election2`) and (iv) deaths (column K of `votes.election2`).
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
#'                      consistent using the iterative proportional fitting algorithm or the Thomsen iterative algorithm.
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
#'  \item{correlations}{ A matrix of order JxK (RxC) with the across units correlations between options for the proportions
#'                       in the transformed scale.}
#'  \item{reference.outputs}{A list with four components: `vjk.averages`, `vjk.units.averages`, `vjk.by.reference`,
#'                           and `vjk.units.by.reference`.  The component `vjk.averages` is a JxKx8 array containing 
#'                           eight global solutions  for the vote transfer matrix, obtained by combining (with different 
#'                           weights) the solutions derived from all combinations of row and column reference options.
#'                           The component `vjk.units.averages` is a JxKxIx8 array containing the corresponding 
#'                           unit-level solutions. The component `vjk.by.reference` is a JxKx(J1`*`K1) array containing 
#'                           all global solutions obtained by considering every possible combination of row and column 
#'                           reference options. The component `vjk.units.by.reference` is a JxKxIx(J1`*`K1) array containing the 
#'                           unit-level solutions associated with `vjk.by.reference`.
#'                           If either `method = "IPF"` or `reference` is not `NULL`, this component is `NULL`.}
#' 
#'  \item{deterministic.bounds}{ A list of two matrices of order JxK (RxC) and two arrays of order JxKxI (RXCxI) containing 
#'                              for each vote transition the lower and upper allowed proportions given the observed aggregates.}
#'
#'  \item{regions.outputs}{ If `regions = NULL`, this is `NULL`. Otherwise, it is a list of ecolRxC outputs, 
#'                          with one element per subset defined by `regions`.}
#'                          
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'  
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
                    regions = NULL,
                    confidence = NULL,
                    B = 0,
                    Yule.aprox = FALSE,
                    ref.combination = "LRCNV",
                    tol = 0.000001,
                    ...){
  
  argg <- c(as.list(environment()), list(...))
  nothing <- tests_inputs_ecolRxC(argg)
  
  if (is.null(regions) | length(unique(regions)) < 2L){
    output <- ecolRxC_r(votes.election1 = votes.election1,
                        votes.election2 = votes.election2,
                        scale = scale,
                        method = method,
                        local = local,
                        census.changes = census.changes,
                        reference = reference,
                        confidence = confidence,
                        B = B,
                        Yule.aprox = Yule.aprox,
                        ref.combination = ref.combination,
                        tol = tol)
  } else {
    inputs <- list("votes.election1" = votes.election1, "votes.election2" = votes.election2,
                   "scale" = scale, "method" = method, "local" = local, "census.changes" = census.changes[1],
                   "regions" = regions, "reference" = reference, "confidence" = confidence, "B" = B, 
                   "Yule.aprox" = Yule.aprox, "ref.combination" = ref.combination, "tol" = tol)
    
    # Estimates across regions
    groups <- unique(regions)
    n.groups <- length(groups)
    regions.outputs <- vector("list", n.groups)
    names(regions.outputs) <- paste("region = ", groups)
    for (gg in 1L:n.groups){
      x.r <- votes.election1[regions == groups[gg], ]
      y.r <- votes.election2[regions == groups[gg], ]
      regions.outputs[[gg]] <- ecolRxC_r(votes.election1 = x.r,
                                        votes.election2 = y.r,
                                        scale = scale,
                                        method = method,
                                        local = local,
                                        census.changes = census.changes,
                                        reference = reference,
                                        confidence = confidence,
                                        B = B,
                                        Yule.aprox = Yule.aprox,
                                        ref.combination = ref.combination,
                                        tol = tol)
    }
    
    # Dealing with differences by regions in new entries and exits
    max.rows <- max(sapply(regions.outputs, function(x) nrow(x$VTM.votes)))
    max.cols <- max(sapply(regions.outputs, function(x) ncol(x$VTM.votes)))
    r.names <- rownames(regions.outputs[[which.max(sapply(regions.outputs, 
                                                          function(x) nrow(x$VTM.votes)))]]$VTM.votes)
    c.names <- colnames(regions.outputs[[which.max(sapply(regions.outputs, 
                                                          function(x) ncol(x$VTM.votes)))]]$VTM.votes)
    
    # Aggregating global estimates across regions
    VTM.votes <- VTM.votes.global <- VTM.crude.global <- matrix(0, max.rows, max.cols)
    dimnames(VTM.votes) <- dimnames(VTM.votes.global) <- dimnames(VTM.crude.global) <-
      list(r.names, c.names)
    VTM.lower <- VTM.upper <- NULL
    if (!is.null(confidence)){
      VTM.lower <- VTM.upper <- matrix(0, max.rows, max.cols)
      dimnames(VTM.lower) <- dimnames(VTM.upper) <- list(r.names, c.names)
    }
    for (gg in 1L:n.groups){
      temp1 <- expand_matrix(regions.outputs[[gg]]$VTM.votes, 
                             max.rows, max.cols, r.names, c.names)
      temp2 <- expand_matrix(regions.outputs[[gg]]$VTM.crude.global, 
                             max.rows, max.cols, r.names, c.names)
      VTM.crude.global <- VTM.crude.global + t(t(temp2)*rowSums(temp1))
      VTM.votes <- VTM.votes + temp1
      VTM.votes.global <- VTM.votes.global + expand_matrix(regions.outputs[[gg]]$VTM.votes.global, 
                                             max.rows, max.cols, r.names, c.names)
      if (!is.null(confidence)){
        temp2 <- expand_matrix(regions.outputs[[gg]]$VTM.lower, 
                               max.rows, max.cols, r.names, c.names)
        VTM.lower <- VTM.lower + sweep(temp2, 1, rowSums(temp1), "*")
        temp2 <- expand_matrix(regions.outputs[[gg]]$VTM.upper, 
                               max.rows, max.cols, r.names, c.names)
        VTM.upper <- VTM.upper + sweep(temp2, 1, rowSums(temp1), "*")
      }
    }
    VTM <- VTM.votes/rowSums(VTM.votes)
    VTM.global <- VTM.votes.global/rowSums(VTM.votes.global)
    VTM.crude.global <- VTM.crude.global/rowSums(VTM.votes)
    if (!is.null(confidence)){
      VTM.lower <- pmin(VTM, VTM.lower/rowSums(VTM.votes))
      VTM.upper <- pmax(VTM, VTM.upper/rowSums(VTM.votes))
    }

    # Combining local estimates by regions
    VTM.lower.units <- VTM.upper.units <- reference.outputs <- NULL
    if( local | method == "Thomsen"){
      VTM.votes.units <- VTM.crude.units <- VTM.units <- 
           array(0, c(max.rows, max.cols, nrow(votes.election1)))
      dimnames(VTM.votes.units) <- dimnames(VTM.crude.units) <- dimnames(VTM.units) <- 
           list(r.names, c.names, rownames(votes.election1))
      if (!is.null(confidence)){
        VTM.lower.units <- VTM.upper.units <- 
          array(0, c(max.rows, max.cols, nrow(votes.election1)))
        dimnames(VTM.lower.units) <- dimnames(VTM.lower.units) <-  
          list(r.names, c.names, rownames(votes.election1))
      }
      for (gg in 1L:n.groups){
        temp1 <- expand_array(regions.outputs[[gg]]$VTM.votes.units, max.rows, max.cols)
        temp2 <- expand_array(regions.outputs[[gg]]$VTM.crude.units, max.rows, max.cols)
        temp3 <- expand_array(regions.outputs[[gg]]$VTM.units, max.rows, max.cols)
        dimnames(temp1) <- dimnames(temp2) <- list(r.names, c.names, 
                                                   rownames(votes.election1)[regions == groups[gg]])
        VTM.votes.units[, , regions == groups[gg]] <- temp1
        VTM.crude.units[, , regions == groups[gg]] <- temp2
        VTM.units[, , regions == groups[gg]] <- temp3
        if (!is.null(confidence)){
          temp1 <- expand_array(regions.outputs[[gg]]$VTM.lower.units, max.rows, max.cols)
          temp2 <- expand_array(regions.outputs[[gg]]$VTM.upper.units, max.rows, max.cols)
          dimnames(temp1) <- dimnames(temp2) <- list(r.names, c.names, 
                                                     rownames(votes.election1)[regions == groups[gg]])
          VTM.lower.units[, , regions == groups[gg]] <- temp1
          VTM.upper.units[, , regions == groups[gg]] <- temp2
        }
      }
      # Combining reference outputs across regions
      if (is.null(reference) & method == "Thomsen"){
        n.method <- dimnames(regions.outputs[[1L]]$reference.outputs$vjk.averages)[[3]]
        n.cells <- dimnames(regions.outputs[[1L]]$reference.outputs$vjk.by.reference)[[3]]
        vjk.averages <- array(0, c(max.rows, max.cols, 8L))
        vjk.units.averages <- array(0, c(max.rows, max.cols, nrow(votes.election1), 8L))
        vjk.by.reference <- array(0, c(max.rows, max.cols, ncol(votes.election1)*ncol(votes.election2)))
        vjk.units.by.reference <- array(0, c(max.rows, max.cols, nrow(votes.election1),
                                             ncol(votes.election1)*ncol(votes.election2)))
        dimnames(vjk.averages) <- list(r.names, c.names, n.method)
        dimnames(vjk.units.averages) <- list(r.names, c.names, rownames(votes.election1), n.method)
        dimnames(vjk.by.reference) <- list(r.names, c.names, n.cells)
        dimnames(vjk.units.by.reference) <- list(r.names, c.names, rownames(votes.election1), n.cells) 
        for (gg in 1L:n.groups){
          temp1 <- expand_array(regions.outputs[[gg]]$reference.outputs$vjk.averages, 
                                max.rows, max.cols)
          temp2 <- expand_array(regions.outputs[[gg]]$reference.outputs$vjk.by.reference, 
                                max.rows, max.cols)
          dimnames(temp1) <- list(r.names, c.names, n.method)
          dimnames(temp2) <- list(r.names, c.names, n.cells)
          vjk.averages <- vjk.averages + temp1
          vjk.by.reference <- vjk.by.reference + temp2
          temp <- regions.outputs[[gg]]$reference.outputs$vjk.units.averages
          vjk.units.averages[1L:dim(temp)[1L], 1L:dim(temp)[2L], regions == groups[gg], ] <- temp
          temp <- regions.outputs[[gg]]$reference.outputs$vjk.units.by.reference
          vjk.units.by.reference[1L:dim(temp)[1L], 1L:dim(temp)[2L], regions == groups[gg], ] <- temp
        }
        reference.outputs <- list("vjk.averages" = vjk.averages, "vjk.units.averages" = vjk.units.averages,
                                  "vjk.by.reference" = vjk.by.reference,
                                  "vjk.units.by.reference" = vjk.units.by.reference)
      }
    }
    # Correlations
    correlations <- matrix(NA, max.rows, max.cols)
    dimnames(correlations) <- list(r.names, c.names)
    if (scale == "logit"){
      ecol2x2 <- Thomsen_logit_2x2
    } else {
      ecol2x2 <- Thomsen_probit_2x2
    }
    net <- compute_net_voters(x0 = votes.election1, y0 = votes.election2, 
                              scenario = census.changes[1L])
    x <- net$x
    y <- net$y
    if (census.changes[1L] == "adjust"){
      x <- votes.election1 * (rowSums(votes.election2)/rowSums(votes.election1))
      y <- votes.election2
    }
    for (j in 1L:ncol(x)){
      for (k in 1L:ncol(y)){
        mt <- ecol2x2(x[, j], rowSums(x), y[, k], rowSums(y),
                      confidence = NULL, Yule.aprox = Yule.aprox)
        correlations[j, k] <- mt$cor
      }
    }
    # Deterministic bounds
    det.bounds <- bounds_compound(origin = x, destination = y, scenario = census.changes[1L],
                                  J0 = ncol(votes.election1), K0 = ncol(votes.election2))
    # Generating output
    output <- list("VTM" = VTM, "VTM.votes" = VTM.votes, "VTM.global" = VTM, "VTM.votes.global" = VTM.votes,
                   "VTM.lower" = VTM.lower, "VTM.upper" = VTM.upper, "VTM.crude.global" = VTM.crude.global,
                   "VTM.units" = VTM.units, "VTM.votes.units" = VTM.votes.units, 
                   "VTM.lower.units" = VTM.lower.units, "VTM.upper.units" = VTM.upper.units, 
                   "VTM.crude.units" = VTM.crude.units, "correlations" = correlations,
                   "reference.outputs" = reference.outputs, "deterministic.bounds" = det.bounds, 
                   "regions.outputs" = regions.outputs, "inputs" = inputs)
  }
  
  class(output) <- "ecolRxC"
  return(output)
}




ecolRxC_r <- function(votes.election1,
                      votes.election2,
                      scale = "probit",
                      method = "Thomsen",
                      local = TRUE,
                      census.changes = c("adjust", "raw", "regular", "ordinary", "enriched",
                                         "simultaneous", "semifull", "full", "gold"),
                      reference = NULL,
                      confidence = NULL,
                      B = 0,
                      Yule.aprox = FALSE,
                      ref.combination = "LRCNV",
                      tol = 0.000001,
                      ...){

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
                              ref.combination = ref.combination,
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
                          ref.combination = "ACVR",
                          tol = 0.000001,
                          ...){

  # argg <- c(as.list(environment()), list(...))
  x0 <- as.matrix(votes.election1)
  y0 <- as.matrix(votes.election2)


  # inputs
  inputs <- list("votes.election1" = votes.election1, "votes.election2" = votes.election2,
                 "scale" = scale, method = "IPF", "local" = FALSE, "census.changes" = census.changes[1],
                 "confidence" = confidence, "B" = B, "Yule.aprox" = Yule.aprox, "ref.combination" = ref.combination,
                 "tol" = tol)

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
                   "reference.outputs" = NULL, "inputs" = inputs)
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
  # deterministic bounds
  det.bounds <- bounds_compound(origin = x, destination = y, scenario = scenario,
                                J0 = J0, K0 = K0)
  
  # Adjustment of initial solutions for congruence
  VTM.votes <- IPF(VTM.crude, vector.columna, vector.fila, precision = tol)
  VTM <- VTM.votes/rowSums(VTM.votes)

  # Uncertainty
  VTM.lower <- VTM.upper <- NULL
  if(!is.null(confidence)){
    if (B > 0.5){
      muestra <- extract_sample(TM.low = VTM.crude.l, TM.upp = VTM.crude.u, B = B)
      intervalos <- interval_transitions(muestra = muestra,
                                         vector.fila = vector.fila,
                                         vector.columna = vector.columna,
                                         tol = tol,
                                         confidence = confidence)
      VTM.lower <- intervalos$TM.low
      VTM.upper <- intervalos$TM.upp
    } else {
      intervalos <- interval_transitions_conservative(TM.low.c = VTM.crude.l, 
                                                      TM.upp.c = VTM.crude.u,
                                                      vector.fila = vector.fila,
                                                      vector.columna = vector.columna,
                                                      tol = tol)
      VTM.lower <- intervalos$TM.low
      VTM.upper <- intervalos$TM.upp
    }
    VTM.lower <- pmax(VTM.lower, det.bounds$lower)
    VTM.upper <- pmin(VTM.upper, det.bounds$upper)
    
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
                 "reference.outputs" = NULL, "deterministic.bounds" = det.bounds, "inputs" = inputs)
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
                          ref.combination = "ACVR",
                          tol = 0.000001,
                          ...){

  # argg <- c(as.list(environment()), list(...))
  x0 <- as.matrix(votes.election1)
  y0 <- as.matrix(votes.election2)

  # inputs
  inputs <- list("votes.election1" = votes.election1, "votes.election2" = votes.election2,
                 "scale" = scale, method = "IPF", "local" = TRUE, "census.changes" = census.changes[1],
                 "confidence" = confidence, "B" = B, "Yule.aprox" = Yule.aprox, "ref.combination" = ref.combination,
                 "tol" = tol)

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

  # deterministic bounds
  det.bounds <- bounds_compound(origin = x, destination = y, scenario = scenario,
                                J0 = J0, K0 = K0)
  
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
                   "reference.outputs" = NULL, "deterministic.bounds" = det.bounds,
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
    if (B > 0.5){
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
    } else {
      for (ii in 1L:I){
        intervalos <- interval_transfers_conservative(TM.low.c = VTM.crude.l.local[, , ii], 
                                                      TM.upp.c = VTM.crude.u.local[, , ii],
                                                      vector.fila = x[ii, ],
                                                      vector.columna = y[ii, ],
                                                      tol = tol)
        VTM.crude.l.local[, , ii] <- intervalos$TM.low/rowSums(VTM.votes.units[, , ii])
        VTM.crude.u.local[, , ii] <- intervalos$TM.upp/rowSums(VTM.votes.units[, , ii])
        VTM.crude.l.local[x[ii, ] == 0L, , ii] <- 0L
        VTM.crude.u.local[x[ii, ] == 0L, , ii] <- 0L
        VTM.lower <- VTM.lower + intervalos$TM.low
        VTM.upper <- VTM.upper + intervalos$TM.upp
     }
  }
    VTM.lower <- VTM.lower/rowSums(VTM.votes)
    VTM.upper <- VTM.upper/rowSums(VTM.votes)
    VTM.lower <- pmax(VTM.lower, det.bounds$lower)
    VTM.upper <- pmin(VTM.upper, det.bounds$upper)
    VTM.crude.l.local <- pmax(VTM.crude.l.local, det.bounds$lower.units)
    VTM.crude.u.local <- pmin(VTM.crude.u.local, det.bounds$upper.units)
    
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
                 "reference.outputs" = NULL, "deterministic.bounds" = det.bounds, 
                 "inputs" = inputs)

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
                            ref.combination = "LRCNV",
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
                 "B" = B, "Yule.aprox" = Yule.aprox, "ref.combination" = ref.combination, "tol" = tol)

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

  # deterministic bounds
  det.bounds <- bounds_compound(origin = x, destination = y, scenario = scenario,
                                J0 = J0, K0 = K0)
  
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
                   "reference.outputs" = NULL, "deterministic.bounds" = det.bounds, 
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
    # Solution combining reference solutions
    #weights <- abs(as.vector(t(correlations[1L:J0, 1L:K0])))
    #weights <- weights/sum(weights)
    #W <- array(rep(weights, each = I*J*K), dim(vjk.units.multi))
    #vjk.units <- apply(vjk.units.multi * W, c(1L, 2L, 3L), sum)
    vjk.units.averages <- units_by_weights(vjk.units.multi = vjk.units.multi,
                                           VTM.crude = VTM.crude, 
                                           pjk.crude.local = pjk.crude.local,
                                           x0 = x0, y0 = y0, x = x, y= y,
                                           tol = tol, correlations = correlations)
    vjk.units <- vjk.units.averages[, , , dimnames(vjk.units.averages)[[4L]] == ref.combination]
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
    reference.outputs <- reference_outputs(vjk.units.multi = vjk.units.multi, 
                                           vjk.units.averages = vjk.units.averages,
                                           VTM.crude = VTM.crude,
                                           pjk.crude.local = pjk.crude.local, 
                                           x0 = x0, y0 = y0,
                                           correlations = correlations)
  }

  # Uncertainty
  VTM.upper <- VTM.lower <- matrix(NA, J, K)
  if(!is.null(confidence)){
    if (B > 0.5){
      muestra_vjk_local <- array(NA, c(J, K, I, B))
      for (bb in 1L:B){
        muestra_vjk_local[, , , bb] <- extract_muestra_vjk(pjk.low = pjk.crude.l.local,
                                                           pjk.upp = pjk.crude.u.local,
                                                           Yule.aprox = Yule.aprox,
                                                           reference = reference,
                                                           scale = scale, x = x, y = y,
                                                           J0 = J0, K0 = K0, tol = tol,
                                                           correlations = correlations,
                                                           ref.combination = ref.combination,
                                                           VTM.crude = VTM.crude)
      }
      # Confidence Intervals
      int.conf <- intervals_Thomsen(muestra_vjk_local = muestra_vjk_local, vjk.units = vjk.units,
                                    x = x, confidence = confidence)
      pjk.crude.l.local <- int.conf$VTM.l.local
      pjk.crude.u.local <- int.conf$VTM.u.local
      VTM.lower <- int.conf$VTM.lower
      VTM.upper <- int.conf$VTM.upper
    } else { 
      selec_vjk_local_min <- selec_vjk_local_max <- array(NA, c(J, K, I, J*K))
      for (bb in 1L:(J*K)){
        selec_vjk_local_min[, , , bb] <- extract_selec_vjk_min(pjk.low = pjk.crude.l.local,
                                                               pjk.upp = pjk.crude.u.local,
                                                               Yule.aprox = Yule.aprox,
                                                               reference = reference,
                                                               scale = scale, x = x, y = y,
                                                               J0 = J0, K0 = K0, tol = tol,
                                                               correlations = correlations,
                                                               bb = bb,
                                                               ref.combination = ref.combination,
                                                               VTM.crude = VTM.crude)
        
        selec_vjk_local_max[, , , bb] <- extract_selec_vjk_max(pjk.low = pjk.crude.l.local,
                                                               pjk.upp = pjk.crude.u.local,
                                                               Yule.aprox = Yule.aprox,
                                                               reference = reference,
                                                               scale = scale, x = x, y = y,
                                                               J0 = J0, K0 = K0, tol = tol,
                                                               correlations = correlations,
                                                               bb = bb,
                                                               ref.combination = ref.combination,
                                                               VTM.crude = VTM.crude)
        
      }
      # Confidence Intervals
      int.conf <- intervals_Thomsen_conservative(selec_vjk_local_min, selec_vjk_local_max, 
                                                 vjk.units, x)
        
      pjk.crude.l.local <- int.conf$VTM.l.local
      pjk.crude.u.local <- int.conf$VTM.u.local
      VTM.lower <- int.conf$VTM.lower
      VTM.upper <- int.conf$VTM.upper
    }
    
    VTM.lower <- pmax(VTM.lower, det.bounds$lower)
    VTM.upper <- pmin(VTM.upper, det.bounds$upper)
    pjk.crude.l.local <- pmax(pjk.crude.l.local, det.bounds$lower.units)
    pjk.crude.u.local <- pmin(pjk.crude.u.local, det.bounds$upper.units)
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
                 "reference.outputs" = reference.outputs, "deterministic.bounds" = det.bounds,
                 "inputs" = inputs)

  return(output)
}


