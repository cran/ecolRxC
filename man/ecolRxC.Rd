% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecolRxC.R
\name{ecolRxC}
\alias{ecolRxC}
\title{Ecological Inference of RxC Tables by Latent Structure Approaches}
\usage{
ecolRxC(
  votes.election1,
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
  tol = 1e-06,
  ...
)
}
\arguments{
\item{votes.election1}{data.frame (or matrix) of order IxJ1 with the votes gained by
(or the counts corresponding to) the J1 (social classes) political options competing
(available) on election 1 (or origin) in the I units considered.}

\item{votes.election2}{data.frame (or matrix) of order IxK2 with the votes gained by
(or the counts corresponding to) the K2 political options competing
(available) on election 2 (or destination) in the I (territorial) units considered.}

\item{scale}{A character string indicating the type of transformation to be applied to the vote
fractions for applying ecological inference. Only \code{logit} and \code{probit} are allowed.
Default, \code{probit}.}

\item{method}{A character string indicating the algorithm to be used for adjusting (making congruent with
the observed margins) the initial crude fractions attained in a 2x2 fashion.
Only \code{Thomsen} (see sec. 4.3 in Thomsen, 1987) and \code{IPF} (iterative proportional fitting,
also known as raking) are allowed. This argument has no effect in the 2x2 case.
Default, \code{Thomsen}.}

\item{local}{A TRUE/FALSE argument indicating whether local solutions (solutions for each polling unit)
must be computed. In that case, the global solution is attained as composition/aggregation
of local solutions. When \code{method = "Thomsen"} local solutions are always computed.
Default \code{TRUE}.}

\item{census.changes}{A character string informing about the level of information available
in \code{votes.election1} and \code{votes.election2} regarding new entries
and exits of the election censuses between the two elections or
indicating how their sum discrepancies should be handled.
This argument allows the eight options discussed in Pavia (2022)
as well as an adjusting option. This argument admits nine values: \code{adjust},
\code{raw}, \code{regular}, \code{ordinary}, \code{simultaneous}, \code{enriched}, \code{semifull},
\code{full} and \code{gold}. See \strong{Details}. Default, \code{adjust}.}

\item{reference}{A vector of two components indicating (parties) options in election 1
and 2, respectively, to be used as reference with \code{method = "Thomsen"}. This has not effect
with \code{method = "IPF"}. The references can be indicated by name or by position.
If \code{reference = NULL}, the final solution is constructed as a weighted average
of all the congruent solutions attained after considering as references all
combinations of options. Default \code{NULL}.}

\item{confidence}{A number between 0 and 1 to be used as level of confidence for the
confidence intervals of the transition rates. By default \code{NULL}.
If \code{confidence = NULL}, confidence intervals are not computed.}

\item{B}{An integer indicating the number of samples to be drawn from each crude estimated
confidence interval for estimating final confidence intervals when either
R (J) or C (K) is higher than two. This is not relevant for the 2x2 case.
It can take a while to compute confidence intervals, mainly when \code{method = Thomsen}.
In general computation burden grows with \code{B}. Default, \code{500}.}

\item{Yule.aprox}{\code{TRUE}/\code{FALSE} argument indicating if either Thomsen (1987)'s formula (3.44),
based on a binormal, or Thomsen (1987)'s formula (3.46), based on Yule's
approximation of tetrachoric correlation, should be use to estimate
cross-proportions. Default \code{FALSE}, formula (3.44).}

\item{tol}{A number indicating the level of precision to be used to stop the
adjustment of initial/crude count estimates reached using a 2x2 approach in a
general RxC case. This is not relevant for the 2x2 case. Default, \code{0.000001}.}

\item{...}{Other arguments to be passed to the function. Not currently used.}
}
\value{
A list with the following components
\item{VTM}{ A matrix of order JxK (RxC) with the estimated proportions of the row-standardized vote transitions from election 1 to election 2.
In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios, this matrix includes the row and the column corresponding to net entries
and net exits (when they are present). When \code{local = TRUE} (default), this matrix is obtained as
composition of the local solutions.}
\item{VTM.votes}{ A matrix of order JxK (RxC) with the estimated vote transfers from election 1 to election 2.
In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios, this matrix includes the row and the column corresponding to net entries
and net exits (when they are present). When \code{local = TRUE} (default), this matrix is obtained as
aggregation of the local solutions.}
\item{VTM.global}{ A matrix of order JxK (RxC) with the estimated proportions of the row-standardized vote transitions from election 1 to election 2,
attained directly from the global (whole electoral space) proportions. When \code{local = FALSE}. \code{VTM} and \code{VTM.global} coincide.
In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios, this matrix includes the row and the column corresponding to net entries
and net exits (when they are present).}
\item{VTM.votes.global}{ A matrix of order JxK (RxC) with the estimated vote transfers from election 1 to election 2,
attained directly from the global proportions. When \code{local = FALSE}, \code{VTM.votes} and \code{VTM.votes.global} coincide.
In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios, this matrix includes the row and the column corresponding to net entries
and net exits (when they are present).}
\item{VTM.lower}{ A matrix of order JxK (RxC) with the estimated lower limits of the confidence intervals for
the proportions of the row-standardized vote transitions from election 1 to election 2.
In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios, this matrix includes the row and the column corresponding to net entries
and net exits (when they are present). When \code{confidence = NULL} this is a \code{NULL} object.}
\item{VTM.upper}{ A matrix of order JxK (RxC) with the estimated upper limits of the confidence intervals for
the proportions of the row-standardized vote transitions from election 1 to election 2.
In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios, this matrix includes the row and the column corresponding to net entries
and net exits (when they are present). When \code{confidence = NULL} this is a \code{NULL} object.}
\item{VTM.crude.global}{ A matrix of order JxK (RxC) with the  (inconsistent) crude estimated proportions for the row-standardized
vote transitions from election 1 to election 2 in the whole space attained in a 2x2 fashion before making them
consistent using the iterative proportional fitting algorithm or the Thomsen iteratuve algortihm.
In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios, this matrix includes the row and the
column corresponding to net entries and net exits (when they are present).}
\item{VTM.units}{ An array of order JxKxI (RxCxI) with the estimated proportions of the row-standardized vote transitions from election 1 to election 2
attained for each unit. When \code{local = FALSE}, this is a \code{NULL} object.
In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios, each unit matrix includes the row and the column corresponding to net entries
and net exits (when they are present).}
\item{VTM.votes.units}{ An array of order JxKxI (RxCxI) with the estimated transfer of votes from election 1 to election 2
attained for each unit. When \code{local = FALSE}, this is a \code{NULL} object.
In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios, each unit matrix includes the row and the column corresponding to net entries
and net exits (when they are present).}
\item{VTM.lower.units}{ An array of order JxKxI (RxCxI) with the estimated lower limits of the confidence intervals for
the proportions of the row-standardized vote transitions from election 1 to election 2 corresponding to each unit.
When either \code{local = FALSE} or \code{confidence = NULL}, this is a \code{NULL} object.
In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios, each unit matrix includes the row and the column corresponding to net entries
and net exits (when they are present).}
\item{VTM.upper.units}{ An array of order JxKxI (RxCxI) with the estimated upper limits of the confidence intervals for
the proportions of the row-standardized vote transitions from election 1 to election 2 corresponding to each unit.
When either \code{local = FALSE} or \code{confidence = NULL}, this is a \code{NULL} object.
In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios, each unit matrix includes the row and the column corresponding to net entries
and net exits (when they are present).}
\item{VTM.crude.units}{ An array of order JxKxI (RxCxI) with the (inconsistent) crude estimated proportions of the row-standardized vote transitions from election 1 to election 2
attained for each unit in a 2x2 fashion before making them consistent using the iterative proportional fitting algorithm or the Thomsen iterative algorithm.
When \code{local = FALSE}, this is a NULL object. In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios, each unit matrix includes the row and the column
corresponding to net entries and net exits (when they are present).}
\item{correlations}{ A matrix of order JxK (Rxc) with the across units correlations between options for the proportions
in the transformed scale.}
\item{reference.outputs}{ A list with three components: \code{vjk.averages}, \code{vjk.by.reference} and \code{vjk.units.by.reference}.
The first component \code{vjk.averages} is a JxKx8 array with eight different global solutions
of transfer matrix of votes attained after combining with different weights each of the
solutions obtained using the different combinations of a row and a column option as reference.
The second component \code{vjk.by.reference} is a JxKx(J1·K1) array with the J1K1 different
global solutions of transfer matrix of votes attained after choosing as reference all the
possible combination of a row and a column option. The third component \code{vjk.units.by.reference}
is a JxKxIx(J1·K1) array with the local solutions linked to \code{vjk.by.reference}.
When either \code{method = "IPF"} or \code{reference} is not \code{NULL}, this is a \code{NULL} object.}

\item{iter}{ A vector of either length 1 (when \code{reference} is different of \code{NULL}) or J1·K1 with
the number of iterations needed by the Thomsen algorithm to reach convergence for
each reference pair. When \code{method = "Thomsen"} this is a \code{NULL} object.}

\item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
}
\description{
Estimates JxK (RxC) vote transfer matrices (ecological contingency tables) based on Thomsen (1987) and Park (2008) approaches.
}
\details{
Description of the \code{census.changes} argument in more detail.
\itemize{
\item{\code{adjust}: }{The default value. This is the simplest solution for handling discrepancies
between the total number of counts for the first and second elections.
With this value the J1 column-aggregations of the counts
in \code{votes.election1} of the first election are proportionally adjusted to
equal the aggregation of the counts in \code{votes.election2} of the second election.
In this scenario, J is equal to J1 and K equal to K2.}
\item{\code{raw}: }{This argument accounts for a scenario with two elections elapsed at least
some months where only the raw election data recorded in the I (territorial) units,
in which the electoral space under study is divided, are available and net
entries and net exits are approached from the available information.
In this scenario, net exits and net entries are estimated according to
Pavia (2022). When both net entries and exits are no
null, constraint (15) of Pavia (2022) applies: no transfer between entries and
exits are allowed. In this scenario, J could be equal to J1 or J1 + 1 and K equal to
K2 or K2 + 1.}
\item{\code{simultaneous}: }{This is the value to be used in classical ecological inference problems,
such as in ecological studies of social or racial voting, and in scenarios with two simultaneous elections.
In this scenario, the sum by rows of \code{votes.election1} and \code{votes.election2} must coincide.}
\item{\code{regular}: }{This value accounts for a scenario with
two elections elapsed at least some months where (i) the column J1
of \code{votes.election1} corresponds to new (young) electors who have the right
to vote for the first time, (ii) net exits and maybe other additional
net entries are computed according to Pavia (2022). When both net entries and exits
are no null, constraints (13) and (15) of Pavia (2022) apply. In this scenario, J
could be equal to J1 or J1 + 1 and K equal to K2 or K2 + 1.}
\item{\code{ordinary}: }{This value accounts for a scenario
with two elections elapsed at least some months where (i) the column K1
of \code{votes.election2} corresponds to electors who died in the interperiod
election, (ii) net entries and maybe other additional net exits are
computed according to Pavia (2022). When both net entries and net exits are no null,
constraints (14) and (15) of Pavia (2022) apply.
In this scenario, J could be equal to J1 or J1 + 1 and K equal to K2 or K2 + 1.}
\item{\code{enriched}: }{This value accounts for a scenario that somewhat combine \code{regular} and
\code{ordinary} scenarios. We consider two elections elapsed at least some months where
(i) the column J1 of \code{votes.election1} corresponds to new (young) electors
who have the right to vote for the first time, (ii) the column K2 of
\code{votes.election2} corresponds to electors who died in the interperiod
election, (iii) other (net) entries and (net) exits are computed according
to Pavia (2022). When both net entries and net exits are no null, constraints (12) to
(15) of Pavia (2022) apply. In this scenario, J could be equal
to J1 or J1 + 1 and K equal to K2 or K2 + 1.}
\item{\code{semifull}: }{This value accounts for a scenario with two elections elapsed at least some
months, where: (i) the column J1 = J of \code{votes.election1} totals new
electors (young and immigrants) that have the right to vote for the first time in each polling unit and
(ii) the column K2 = K of \code{votes.election2} corresponds to total exits of the census
lists (due to death or emigration). In this scenario, the sum by rows of
\code{votes.election1} and \code{votes.election2} must agree and constraint (15)
of Pavia (2022) apply.}
\item{\code{full}: }{This value accounts for a scenario with two elections elapsed at least some
months, where J = J1, K = K2 and (i) the column J - 1 of \code{votes.election1} totals new (young)
electors that have the right to vote for the first time, (ii) the column J
of \code{votes.election1} measures new immigrants that have the right to vote and
(iii) the column K of \code{votes.election2} corresponds to total exits of the census
lists (due to death or emigration). In this scenario, the sum by rows of
\code{votes.election1} and \code{votes.election2} must agree and constraints (13)
and (15) of Pavia (2022) apply.}
\item{\code{gold}: }{This value accounts for a scenario similar to \code{full}, where J = J1, K = K2 and
total exits are separated out between exits due to emigration
(column K - 1 of \code{votes.election2}) and death (column K of \code{votes.election2}).
In this scenario, the sum by rows of \code{votes.election1} and \code{votes.election2} must agree.
Constraints (12) to (15) of Pavia (2022) apply.}
}
}
\note{
This function somewhere builds on the .ado (STATA) functions written by Won-ho Park, in 2002.
}
\examples{
votes1 <- structure(list(P1 = c(16L, 4L, 13L, 6L, 1L, 16L, 6L, 17L, 48L, 14L),
                         P2 = c(8L, 3L, 0L, 5L, 1L, 4L, 7L, 6L, 28L, 8L),
                         P3 = c(38L, 11L, 11L, 3L, 13L, 39L, 14L, 34L, 280L, 84L),
                         P4 = c(66L, 5L, 18L, 39L, 30L, 57L, 35L, 65L, 180L, 78L),
                         P5 = c(14L, 0L, 5L, 2L, 4L, 21L, 6L, 11L, 54L, 9L),
                         P6 = c(8L, 2L, 5L, 3L, 0L, 7L, 7L, 11L, 45L, 17L),
                         P7 = c(7L, 3L, 5L, 2L, 3L, 17L, 7L, 13L, 40L, 8L)),
                         row.names = c(NA, 10L), class = "data.frame")
votes2 <- structure(list(C1 = c(2L, 1L, 2L, 2L, 0L, 4L, 0L, 4L, 19L, 14L),
                         C2 = c(7L, 3L, 1L, 7L, 2L, 5L, 3L, 10L, 21L, 6L),
                         C3 = c(78L, 7L, 28L, 42L, 28L, 84L, 49L, 85L, 260L, 100L),
                         C4 = c(56L, 14L, 20L, 7L, 19L, 54L, 22L, 50L, 330L, 91L),
                         C5 = c(14L, 3L, 6L, 2L, 3L, 14L, 8L, 8L, 45L, 7L)),
                         row.names = c(NA, 10L), class = "data.frame")
example <- ecolRxC(votes1, votes2, method = "IPF")$VTM

}
\references{
Achen, C.H. (2000). The Thomsen Estimator for Ecological Inference (Unpublished manuscript). University of Michigan.

Park, W.-H. (2008). Ecological Inference and Aggregate Analysis of Elections. PhD Dissertation. University of Michigan.

Pavia, J.M. (2022). Adjustment of initial estimates of voter transition probabilities to guarantee consistency and completeness.

Thomsen, S.R. (1987). Danish Elections 1920-79: a Logit Approach to Ecological Analysis and Inference. Politica, Aarhus, Denmark.
}
\author{
Jose M. Pavia, \email{pavia@uv.es}
}
\concept{latent structure ecological inference functions}
