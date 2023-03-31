Thomsen_logit_2x2 <- function(votes1, total1, votes2, total2, confidence = 0.95,
                              Yule.aprox = FALSE){
# Inputs:
# votes1: matrix/vector nx1 with the votes in the first election of the party of interest
# total1: matrix/vector nx1 with the total votes recorded in the first election
# votes2: matrix/vector nx1 with the votes in the second election of the party of interest
# total2: matrix/vector nx1 with the total votes recorded in the second election
# confidence: level of confidence for the interval estimates. If NULL: not computed.
# Yule.aprox = TRUE/FALSE argument indicating if either formula (3.44) or (3.46) of Thomsen (1987) should be use to estimate cross-proportions. Default FALSE, formula (3.44).
# The variables should be in votes (counts)
# This function estimates the transfer between party1 and party2
# Output: A list with the following elements
# PTM: A 2x2 matrix with the transitions probabilities between votes and total
# PTM_low: A 2x2 matrix with the upper bounds of transitions probabilities between votes and total
# PTM_high: A 2x2 matrix with the lower bounds of transitions probabilities between votes and total
# cor: An estimate of correlation between p(votes1) y p(votes2) in the transformed scale.

#  This function is based on the .m (MATLAB) and .ado (STATA)
#  functions written by Won-ho Park, University of Michigan in July-Aug. 2002
  votes1 <- as.vector(votes1)
  votes2 <- as.vector(votes2)
  total1 <- as.vector(total1)
  total2 <- as.vector(total1)

  I <- length(votes2)
  weight <- total2/sum(total2)
  t1 <- sum(votes1)/sum(total1)
  t2 <- sum(votes2)/sum(total2)

# To avoid taking logs from zero values
  votes1[votes1 == 0L] <- 0.5
  votes2[votes2 == 0L] <- 0.5
  total1[total1 == votes1] <- total1[total1 == votes1] + 0.5
  total2[total2 == votes2] <- total2[total2 == votes2] + 0.5

# Logit transformation of the proportions
  lv1 <- log(votes1/(total1 - votes1))
  lv2 <- log(votes2/(total2 - votes2))
  meanlv1 <- sum(weight * lv1)
  meanlv2 <- sum(weight * lv2)
  varlv1 <- sum(I*weight * (lv1 - meanlv1)^2L)/(I - 1L)
  varlv2 <- sum(I*weight * (lv2 - meanlv2)^2L)/(I - 1L)

# Pearson correlation of logit transformation
  r <- sum(weight * (lv1 - meanlv1) * (lv2 - meanlv2))/sqrt(varlv1*varlv2)

# Integration binormal
  core <- pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r)
  if (Yule.aprox){
    # Tetrachoric formula
    k <- sqrt((1L + 2L*r*t1 + 2L*r*t2 - r)^2L - 8L*r*(1L + r)*t1*t2)
    core <- (1L + 2L*r*t1 + 2L*r*t2 - r - k)/(4L * r)
  }
  p2p <-  as.numeric(core)/t1  # Estimate of vote transfer from party 1 to party 2
  t2p <- (t2 - core)/(1L - t1) # Estimate of vote transfer from total(rest) 1 to party 2
  trans.matrix <- matrix(c(p2p, t2p, 1L - p2p, 1L - t2p), 2L, 2L)
  trans.matrix_low <- trans.matrix_high <- NULL

# Standard deviations
  if (!is.null(confidence)) {
    z <- stats::qnorm(1L - (1L - confidence)/2L)
  # z_hi <- log((1L + r)/(1L - r))/2L + z/sqrt(I - 3L)
  # z_low <- log((1L + r)/(1 - r))/2L - z/sqrt(I - 3L)
  # r_hi <- (exp(2L * z_hi) - 1L)/(exp(2L * z_hi) + 1L)
  # r_low <- (exp(2L * z_low) - 1L)/(exp(2L * z_low) +1L)
    r_hi <- tanh(log((1L + r)/(1L - r))/2L + z/sqrt(I - 3L))
    r_low <- tanh(log((1L + r)/(1L - r))/2L - z/sqrt(I - 3L))

    p2p_hi <- as.numeric(pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r_hi))/t1
    p2p_low <- as.numeric(pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r_low))/t1
    if (Yule.aprox){
      # Tetrachoric formula
      k <- sqrt((1L + 2L*r_hi*t1 + 2L*r_hi*t2 - r_hi)^2L - 8L*r_hi*(1L + r_hi)*t1*t2)
      p2p_hi <- (1L + 2L*r_hi*t1 + 2L*r_hi*t2 - r_hi - k)/(4L * r_hi)/t1
      k <- sqrt((1L + 2L*r_low*t1 + 2L*r_low*t2 - r_low)^2L - 8L*r_low*(1L + r_low)*t1*t2)
      p2p_low <- (1L + 2L*r_low*t1 + 2L*r_low*t2 - r_low - k)/(4L * r_low)/t1
    }
    t2p_hi <- (t2 - p2p_low*t1)/(1L - t1)
    t2p_low <- (t2 - p2p_hi*t1)/(1L - t1)

    trans.matrix_low <- matrix(c(p2p_low, t2p_low, 1L - p2p_hi, 1L - t2p_hi), 2L, 2L)
    trans.matrix_high <- matrix(c(p2p_hi, t2p_hi, 1L - p2p_low, 1L - t2p_low), 2L, 2L)
  }

  output <- list("PTM" = trans.matrix, "PTM_low" = trans.matrix_low,
                 "PTM_high" = trans.matrix_high, "cor" = r)
  return(output)
}

#~~FUNCTIONS FOR CUMULATIVE BIVARIATE NORMAL~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pnorm2d <-function(x, y = x, rho = 0) {
    # pnorm2d: A copy from R package "sn"

    # Description:
    #   Computes bivariate Normal probability function

    # Arguments:
    #   x, y - two numeric values or vectors of the same length at
    #       which the probability will be computed.

    # Value:
    #   returns a numeric vector of probabilities of the same length
    #   as the input vectors

    # FUNCTION:

    # Probaility:
    X <- cbind(x, y)
    ans <- apply(X, 1L, .pnorm2d, rho = rho)
    attr(ans, "control") = c(rho = rho)


    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------

.pnorm2d <- function(X, rho = 0) {
    # .pnorm2d: A copy from R package "sn"

    # Description:
    #   Bivariate Normal probability function

    # Arguments:
    #   x, y - two numeric values at which the probability will
    #   be computed.

    # Value:
    #   returns a numeric vector of probabilities of the same length
    #   as the input vectors

    # FUNCTION:

    # Probability:
    x <- X[1L]
    y <- X[2L]
    if (x == 0L & y == 0L) {
      return(0.25 + asin(rho)/(2L * pi))
    }
    p = 0.5 * (stats::pnorm(x) + stats::pnorm(y))
    if (x == 0L) {
      p = p - 0.25 * sign(y)
    } else {
      if (is.finite(x)) {
        Y = (y - rho * x)/(x * sqrt(1L - rho^2L))
      } else {
        Y = -rho/sqrt(1L - rho^2L)
      }
      p = p - .TOwen(x, Y)
    }
    if (y == 0L) {
      p = p - 0.25 * sign(x)
    } else {
      if (is.finite(y)) {
        X = (x - rho * y)/(y * sqrt(1L - rho^2L))
      } else {
        X = -rho/sqrt(1L - rho^2L)
      }
      p = p - .TOwen(y, X)
    }
    if (is.finite(x) & is.finite(y)) {
      if ((x * y < 0L) | ((x * y == 0L) & (x + y) < 0L)) {
        p = p - 0.5
      }
    }

    # Return Value:
    return(p)
  }


# ------------------------------------------------------------------------------


.TInt <- function(h, a, jmax, cut.point) {
    # .Tint: A copy from R package "sn"

    # Note:
    #   Required by .pnorm2d and .TOwen

    # FUNCTION:

    .fui = function(h, i) (h^(2L * i))/((2L^i) * gamma(i + 1L))
    seriesL = seriesH = NULL
    i = 0L:jmax
    low = (h <= cut.point)
    hL = h[low]
    hH = h[!low]
    L = length(hL)
    if (L > 0L) {
      b = outer(hL, i, .fui)
      cumb = apply(b, 1L, cumsum)
      b1 = exp(-0.5 * hL^2L) * t(cumb)
      matr = matrix(1L, jmax + 1L, L) - t(b1)
      jk = rep(c(1L, -1L), jmax)[1L:(jmax + 1L)]/(2L * i + 1)
      matr = t(matr * jk) %*% a^(2L * i + 1L)
      seriesL = (atan(a) - as.vector(matr))/(2L * pi)
    }
    if (length(hH) > 0L) {
      seriesH = atan(a) * exp(-0.5 * (hH^2L) * a/atan(a)) *
        (1L + 0.00868 * (hH^4L) * a^4L)/(2L * pi)
    }
    series = c(seriesL, seriesH)
    id = c((1L:length(h))[low], (1L:length(h))[!low])
    series[id] = series

    # Return Value:
    series
  }

# ------------------------------------------------------------------------------

.TOwen <- function (h, a, jmax = 50, cut.point = 6) {
    # .TOwen: A copy from R package "sn"

    # Note:
    #   Required by .pnorm2d

    # FUNCTION:

    if (!is.vector(a) | length(a) > 1L)
      stop("a must be a vector of length 1")
    if (!is.vector(h))
      stop("h must be a vector")
    aa = abs(a)
    ah = abs(h)
    if (aa == Inf)
      return(0.5 * stats::pnorm(-ah))
    if (aa == 0L)
      return(rep(0L, length(h)))
    na = is.na(h)
    inf = (ah == Inf)
    ah = replace(ah, (na | inf), 0L)
    if (aa <= 1L) {
      owen = .TInt(ah, aa, jmax, cut.point)
    } else {
      owen = 0.5 * stats::pnorm(ah) + stats::pnorm(aa * ah) * (0.5 - stats::pnorm(ah)) -
        .TInt(aa * ah, (1/aa), jmax, cut.point)
    }
    owen = replace(owen, na, NA)
    owen = replace(owen, inf, 0L)
    ans = return(owen * sign(a))

    # Return Value:
    ans
  }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Thomsen_probit_2x2 <- function(votes1, total1, votes2, total2, confidence = 0.95,
                               Yule.aprox = FALSE){
  # Inputs:
  # votes1: matrix/vector nx1 with the votes in the first election of the party of interest
  # total1: matrix/vector nx1 with the total votes recorded in the first election
  # votes2: matrix/vector nx1 with the votes in the second election of the party of interest
  # total2: matrix/vector nx1 with the total votes recorded in the second election
  # confidence: level of confidence for the interval estimates. If NULL: not computed.
  # Yule.aprox = TRUE/FALSE argument indicating if either formula (3.44) or (3.46) of Thomsen (1987) should be use to estimate cross-proportions. Default FALSE, formula (3.44).
  # The variables should be in votes (counts)
  # This function estimates the transfer between party1 and party2
  # Output: A list with the following elements
  # PTM: A 2x2 matrix with the transitions probabilities between votes and total
  # PTM_low: A 2x2 matrix with the upper bounds of transitions probabilities between votes and total
  # PTM_high: A 2x2 matrix with the lower bounds of transitions probabilities between votes and total
  # cor: An estimate of correlation between p(votes1) y p(votes2) in the transformed scale.

  #  This function is based on the .m (MATLAB) and .ado (STATA)
  #  functions written by Won-ho Park, University of Michigan in July-Aug. 2002
  votes1 <- as.vector(votes1)
  votes2 <- as.vector(votes2)
  total1 <- as.vector(total1)
  total2 <- as.vector(total1)

  I <- length(votes2)
  weight <- total2/sum(total2)
  t1 <- sum(votes1)/sum(total1)
  t2 <- sum(votes2)/sum(total2)

  # To avoid taking probits from one or zero values
  votes1[votes1 == 0L] <- 0.5
  votes2[votes2 == 0L] <- 0.5
  total1[total1 == votes1] <- total1[total1 == votes1] + 0.5
  total2[total2 == votes2] <- total2[total2 == votes2] + 0.5

  # Probit transformation of the proportions
  pv1 <- stats::qnorm(votes1/total1)
  pv2 <- stats::qnorm(votes2/total2)
  meanpv1 <- sum(weight * pv1)
  meanpv2 <- sum(weight * pv2)
  varpv1 <- sum(I*weight * (pv1 - meanpv1)^2L)/(I - 1L)
  varpv2 <- sum(I*weight * (pv2 - meanpv2)^2L)/(I - 1L)

  # Pearson correlation of probit transformations
  r <- sum(weight * (pv1 - meanpv1) * (pv2 - meanpv2))/sqrt(varpv1*varpv2)

  # Integration binormal
  core <- pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r)
  if (Yule.aprox){
    # Tetrachoric formula
    k <- sqrt((1L + 2L*r*t1 + 2L*r*t2 - r)^2L - 8L*r*(1L + r)*t1*t2)
    core <- (1L + 2L*r*t1 + 2L*r*t2 - r - k)/(4L * r)
  }
  p2p <-  as.numeric(core)/t1  # Estimate of vote transfer from party 1 to party 2
  t2p <- (t2 - core)/(1L - t1) # Estimate of vote transfer from total(rest) 1 to party 2
  trans.matrix <- matrix(c(p2p, t2p, 1L - p2p, 1L - t2p), 2L, 2L)
  trans.matrix_low <- trans.matrix_high <- NULL

  # Standard deviations
  if (!is.null(confidence)) {
    z <- stats::qnorm(1L - (1L - confidence)/2L)
    # z_hi <- log((1L + r)/(1L - r))/2L + z/sqrt(I - 3L)
    # z_low <- log((1L + r)/(1 - r))/2L - z/sqrt(I - 3L)
    # r_hi <- (exp(2L * z_hi) - 1L)/(exp(2L * z_hi) + 1L)
    # r_low <- (exp(2L * z_low) - 1L)/(exp(2L * z_low) +1L)
    r_hi <- tanh(log((1L + r)/(1L - r))/2L + z/sqrt(I - 3L))
    r_low <- tanh(log((1L + r)/(1L - r))/2L - z/sqrt(I - 3L))

    p2p_hi <- as.numeric(pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r_hi))/t1
    p2p_low <- as.numeric(pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r_low))/t1
    if (Yule.aprox){
      # Tetrachoric formula
      k <- sqrt((1L + 2L*r_hi*t1 + 2L*r_hi*t2 - r_hi)^2L - 8L*r_hi*(1L + r_hi)*t1*t2)
      p2p_hi <- (1L + 2L*r_hi*t1 + 2L*r_hi*t2 - r_hi - k)/(4L * r_hi)/t1
      k <- sqrt((1L + 2L*r_low*t1 + 2L*r_low*t2 - r_low)^2L - 8L*r_low*(1L + r_low)*t1*t2)
      p2p_low <- (1L + 2L*r_low*t1 + 2L*r_low*t2 - r_low - k)/(4L * r_low)/t1
    }
    t2p_hi <- (t2 - p2p_low*t1)/(1L - t1)
    t2p_low <- (t2 - p2p_hi*t1)/(1L - t1)

    trans.matrix_low <- matrix(c(p2p_low, t2p_low, 1L - p2p_hi, 1L - t2p_hi), 2L, 2L)
    trans.matrix_high <- matrix(c(p2p_hi, t2p_hi, 1L - p2p_low, 1L - t2p_low), 2L, 2L)
  }

  output <- list("PTM" = trans.matrix, "PTM_low" = trans.matrix_low,
                 "PTM_high" = trans.matrix_high, "cor" = r)
  return(output)
}


#-------------------------------------------------------------------------------
# Funcion de ajuste por balanceo
IPF <- function(matriz, vector.columna, vector.fila, precision = 0.000001){
  #
  # IPF(matriz, vector.columna, vector.fila, precision = 0.000001)
  #
  # Funcion para ajustar los coeficientes de una matriz para que verifique las
  # restricciones de agregacion por filas (vector.fila) y columnas (vector.columna)
  # mediante el metodo RAS o the iterative proportional fitting algorithm
  #
  # INPUT:
  #       matriz: matriz mxn inicial (no necesariamente cuadrada) a ajustar/balancear
  #               para que cumpla las restricciones de agregacion por filas y columnas
  #               dadas por los vectores fila y columna
  #       vector.fila: vector de m componentes con lo que debe sumar la matriz
  #                    despues del ajuste
  #       vector.columna: vector de n componentes con la restricciones de lo que deben
  #                       sumar las columas de matriz tras el ajuste
  #       precision: Discrepancia maxima tolerable para alcanzar la convergencia. Por
  #                  defecto 0.000001.
  #
  # OUTPUT:
  #
  #       Salida: matriz mxn balanceada a partir de la matriz inicial cuya suma por filas
  #               coincide con vector.fila y cuya suma por columnas coincide con vector.columna
  #
  # REFERENCIA:
  # Pavia et al (2009): Updating Input-Output Matrices: Assessing
  # Alternatives through Simulation, Journal of Statistical Computation and Simulation.
  #
  # OBSERVACION:
  # El programa no tiene mensajes de aviso por si no se verifica algunas de las
  # restricciones como que sum(vector.fila)=sum(vector.columna)

  nc <- length(vector.columna)
  nf <- length(vector.fila)
  vector.fila1 <- rowSums(matriz)
  R1 <- diag(as.vector(vector.fila)) %*% (diag(as.vector(1/vector.fila1)))
  R1[is.nan(R1)] <- 0L
  R1[is.infinite(R1)] <- 1L
  # R1 <- diag(as.vector(vector.fila)) %*% MASS::ginv(diag(as.vector(vector.fila1)))
  X1 <- R1 %*% matriz
  vector.columna1 <- colSums(X1)
  S1 <- diag(as.vector(vector.columna)) * (diag(as.vector(1/vector.columna1)))
  S1[is.nan(S1)] <- 0L
  S1[is.infinite(S1)] <- 1L
  # S1 <-  diag(as.vector(vector.columna)) * ginv(diag(as.vector(vector.columna1)))
  X2 <- X1 %*% S1
  dif <- max(abs(X2 - matriz))

  while (dif > precision){
    matriz <- X2
    vector.fila1 <- rowSums(matriz)
    R1 <- diag(as.vector(vector.fila)) %*% (diag(as.vector(1/vector.fila1)))
    R1[is.nan(R1)] <- 0L
    R1[is.infinite(R1)] <- 1L
    #   R1 <- diag(as.vector(vector.fila)) %*% ginv(diag(as.vector(vector.fila1)))
    X1 <- R1 %*% matriz
    vector.columna1 <- colSums(X1)
    S1 <- diag(as.vector(vector.columna)) * (diag(as.vector(1/vector.columna1)))
    S1[is.nan(S1)] <- 0L
    S1[is.infinite(S1)] <- 1L
    #   S1 <-  diag(as.vector(vector.columna)) * ginv(diag(as.vector(vector.columna1)))
    X2 <- X1 %*% S1
    dif <- max((abs(X2 - matriz)))
  }
    X2 <- X2*(vector.fila/rowSums(X2))
    X2[is.nan(X2)] <- 0L
    X2[is.infinite(X2)] <- 1L
    X2 <- t(t(X2)*(vector.columna/colSums(X2)))
    X2[is.nan(X2)] <- 0L
    X2[is.infinite(X2)] <- 1L

  return(X2)
}

#-------------------------------------------------------------------------------
# Funcion de ajuste por balanceo 2
IPF2 <- function(matriz, vector.columna, vector.fila, precision = 0.01){

  nc <- length(vector.columna)
  nf <- length(vector.fila)
  vector.fila1 <- rowSums(matriz)
  R1 <- diag(as.vector(vector.fila)) %*% (diag(as.vector(1/vector.fila1)))
  R1[is.nan(R1)] <- 0L
  R1[is.infinite(R1)] <- 1L
  # R1 <- diag(as.vector(vector.fila)) %*% MASS::ginv(diag(as.vector(vector.fila1)))
  X1 <- R1 %*% matriz
  vector.columna1 <- colSums(X1)
  S1 <- diag(as.vector(vector.columna)) * (diag(as.vector(1/vector.columna1)))
  S1[is.nan(S1)] <- 0L
  S1[is.infinite(S1)] <- 1L
  # S1 <-  diag(as.vector(vector.columna)) * ginv(diag(as.vector(vector.columna1)))
  X2 <- X1 %*% S1
  dif.r <- sum(abs(vector.fila - rowSums(X2)))
  dif.c <- sum(abs(vector.columna - colSums(X2)))
  dif <- max(dif.r, dif.c)
  iter <- 0L
  while (dif > precision & iter < 1000L){
    matriz <- X2
    vector.fila1 <- rowSums(matriz)
    R1 <- diag(as.vector(vector.fila)) %*% (diag(as.vector(1/vector.fila1)))
    R1[is.nan(R1)] <- 0L
    R1[is.infinite(R1)] <- 1L
    #   R1 <- diag(as.vector(vector.fila)) %*% ginv(diag(as.vector(vector.fila1)))
    X1 <- R1 %*% matriz
    vector.columna1 <- colSums(X1)
    S1 <- diag(as.vector(vector.columna)) * (diag(as.vector(1/vector.columna1)))
    S1[is.nan(S1)] <- 0L
    S1[is.infinite(S1)] <- 1L
    #   S1 <-  diag(as.vector(vector.columna)) * ginv(diag(as.vector(vector.columna1)))
    X2 <- X1 %*% S1
    dif.r <- sum(abs(vector.fila - rowSums(X2)))
    dif.c <- sum(abs(vector.columna - colSums(X2)))
    dif <- max(dif.r, dif.c)
    iter <- iter + 1L
  }

  return(X2)
}

#-------------------------------------------------------------------------------
# Function that tests if all inputs of ecolRxC are correct and converts
# votes_election1 and votes_election2 into matrix objects
tests_inputs_ecolRxC <- function(argg){

  argg.1 <- as.matrix(argg$votes.election1)
  argg.2 <- as.matrix(argg$votes.election2)

  # Test reference
  if(!is.null(argg$reference)){
    if( !(argg$reference[1] %in% colnames(argg.1) | as.numeric(argg$reference[1]) <= ncol(argg.1)) )
      stop('The value included as first component of "reference" is invalid.')
    if( !(argg$reference[2] %in% colnames(argg.2) | as.numeric(argg$reference[2]) <= ncol(argg.2)) )
      stop('The value included as second component of "reference" is invalid.')
  }

  # Test names 1
  if (!(argg$census.changes[1L] %in% c("adjust", "raw", "regular", "ordinary", "semifull", "enriched",
                                            "simultaneous", "full", "fullreverse", "gold")))
    stop('The value set for argument "census.changes" is not allowed. The only allowed strings for "census.changes" are "adjust", "raw", "simultaneous", "regular", "ordinary", "enriched", "semifull", "full" and "gold".')

  # Test names 2
  if (!(argg$scale[1L] %in% c("logit", "probit")))
    stop('The value set for argument "scale" is not allowed. The only allowed strings for "scale" are "logit" and "probit".')

  # Test method
  if (!(argg$method[1L] %in% c("Thomsen", "IPF")))
    stop('The value set for argument "method" is not allowed. The only allowed strings for "method" are "Thomsen" and "IPF".')


  # Test logical
  if (!is.logical(argg$local))
    stop('The value set for argument "local" is not allowed. The only allowed values are the logical values "TRUE" or "FALSE".')

  # Tests numeric data
  x <- as.matrix(argg$votes.election1)
  y <- as.matrix(argg$votes.election2)
  if (ncol(x) < 2L)
    stop('At least to options/groups/parties/columns are required in "votes.election1".')
  if (ncol(y) < 2L)
    stop('At least to options/groups/parties/columns are required in "votes.election2".')
  if (nrow(x) != nrow(y))
    stop('The number of (spatial) units is different in origin and destination.')
  if (argg$census.changes[1L] %in% c("semifull", "simultaneous", "full", "fullreverse", "gold")){
    if (!identical(round(rowSums(x)), round(rowSums(y)))){
      texto <- paste0('The number of voters (electors) in Election 1 and ',
                      'Election 2 differ (in at least a territorial unit). This is not ',
                      'allowed in a \"', argg$census.changes[1L], '\" scenario.')
      stop(texto)
    }
  }
  if (min(x,y) < 0L) stop('Negative counts are not allowed in arguments "votes.election1" and "votes.election2".')

  # Test confidence
  if(!is.null(argg$confidence)){
    if(argg$confidence > 1L | argg$confidence < 0L)
      stop('Only values in the interval (0, 1) are allowed for the "confidence" argument.')
  }

  # Test B
  dif <- argg$B - round(argg$B)
  if(dif != 0L | argg$B < 0L)
      stop('Only positive integer values are allowed for the "B" argument.')

  # Test tol
  if(argg$tol < 0L)
  stop('Only positive values are allowed for the "tol" argument.')

  return(list("x" = x, "y" = y))
}


#-------------------------------------------------------------------------------
# Function to compute net entries and exits and to generate the origin and destination matrices,
# including if it is the case net entries and exits
compute_net_voters <- function(x0, y0, scenario){
  NET_ENTRIES <- NET_EXITS <- rep(0L, nrow(x0))
  x <- x0
  y <- y0

  if (scenario != "adjust"){
    # Estimation of net entries/exits in the election census
    d <- rowSums(y) - rowSums(x)
    if (any(d != 0L)) {
      NET_ENTRIES[d > 0L] <- d[d > 0L]
      NET_EXITS[d < 0L] <- -d[d < 0L]
    }

    # Net entries and exits
    if (sum(NET_ENTRIES) > 0L){
      x <- cbind(x, NET_ENTRIES)
    }
    if (sum(NET_EXITS) > 0L){
      y <- cbind(y0, NET_EXITS)
    }
  }

  return(list("x" = x, "y" = y))
}

#-------------------------------------------------------------------------------

# Function to add the constraints given by the scenario

constraints_scenario <- function(VTM.crude, VTM.crude.l, VTM.crude.u, scenario, J0, K0){
  J <- nrow(VTM.crude)
  K <- ncol(VTM.crude)
  if((scenario == "raw" & J > J0 & K > K0) |
     (scenario == "regular" & J == J0 & K > K0) |
     (scenario == "ordinary" & J > J0 & K == K0) |
     (scenario == "enriched" & J == J0 & K == K0) |
     (scenario == "semifull")) {
        VTM.crude[J, K] <- VTM.crude.l[J, K] <- VTM.crude.u[J, K] <- 0L
  } else if ((scenario == "regular" & J > J0 & K > K0) |
             (scenario == "enriched" & J > J0 & K == K0) |
             (scenario == "full")) {
    VTM.crude[J, K] <- VTM.crude.l[J, K] <- VTM.crude.u[J, K] <- 0L
    VTM.crude[J - 1L, K] <- VTM.crude.l[J - 1L, K] <- VTM.crude.u[J - 1L, K] <- 0L
  } else if ((scenario == "ordinary" & J > J0 & K > K0) |
             (scenario == "enriched" & J == J0 & K > K0)) {
    VTM.crude[J, K] <- VTM.crude.l[J, K] <- VTM.crude.u[J, K] <- 0L
    VTM.crude[J, K - 1L] <- VTM.crude.l[J, K - 1L] <- VTM.crude.u[J, K -1L] <- 0L
  } else if ((scenario == "enriched" & J > J0 & K > K0) |
             (scenario == "gold")) {
    VTM.crude[J, K] <- VTM.crude.l[J, K] <- VTM.crude.u[J, K] <- 0L
    VTM.crude[J, K - 1L] <- VTM.crude.l[J, K - 1L] <- VTM.crude.u[J, K -1L] <- 0L
    VTM.crude[J - 1L, K] <- VTM.crude.l[J - 1L, K] <- VTM.crude.u[J - 1L, K] <- 0L
    VTM.crude[J - 1L, K - 1L] <- VTM.crude.l[J - 1L, K - 1L] <- VTM.crude.u[J - 1L, K - 1L] <- 0L
  }
  output <- list("VTM.crude" = VTM.crude, "VTM.crude.l" = VTM.crude.l, "VTM.crude.u" = VTM.crude.u)
  return(output)
}

#-------------------------------------------------------------------------------

# Function to add the constraints given by the scenario at the unit level

constraints_scenario_local <- function(VTM.crude.local, VTM.crude.l.local, VTM.crude.u.local, scenario, J0, K0){
  J <- dim(VTM.crude.local)[1L]
  K <- dim(VTM.crude.local)[2L]
  if((scenario == "raw" & J > J0 & K > K0) |
     (scenario == "regular" & J == J0 & K > K0) |
     (scenario == "ordinary" & J > J0 & K == K0) |
     (scenario == "enriched" & J == J0 & K == K0) |
     (scenario == "semifull")) {
    VTM.crude.local[J, K, ] <- VTM.crude.l.local[J, K, ] <- VTM.crude.u.local[J, K, ] <- 0L
  } else if ((scenario == "regular" & J > J0 & K > K0) |
             (scenario == "enriched" & J > J0 & K == K0) |
             (scenario == "full")) {
    VTM.crude.local[J, K, ] <- VTM.crude.l.local[J, K, ] <- VTM.crude.u.local[J, K, ] <- 0L
    VTM.crude.local[J - 1L, K, ] <- VTM.crude.l.local[J - 1L, K, ] <- VTM.crude.u.local[J - 1L, K, ] <- 0L
  } else if ((scenario == "ordinary" & J > J0 & K > K0) |
             (scenario == "enriched" & J == J0 & K > K0)) {
    VTM.crude.local[J, K, ] <- VTM.crude.l.local[J, K, ] <- VTM.crude.u.local[J, K, ] <- 0L
    VTM.crude.local[J, K - 1L, ] <- VTM.crude.l.local[J, K - 1L, ] <- VTM.crude.u.local[J, K -1L, ] <- 0L
  } else if ((scenario == "enriched" & J > J0 & K > K0) |
             (scenario == "gold")) {
    VTM.crude.local[J, K, ] <- VTM.crude.l.local[J, K, ] <- VTM.crude.u.local[J, K, ] <- 0L
    VTM.crude.local[J, K - 1L, ] <- VTM.crude.l.local[J, K - 1L, ] <- VTM.crude.u.local[J, K -1L, ] <- 0L
    VTM.crude.local[J - 1L, K, ] <- VTM.crude.l.local[J - 1L, K, ] <- VTM.crude.u.local[J - 1L, K, ] <- 0L
    VTM.crude.local[J - 1L, K - 1L, ] <- VTM.crude.l.local[J - 1L, K - 1L, ] <- VTM.crude.u.local[J - 1L, K - 1L, ] <- 0L
  }
  output <- list("VTM.crude.local" = VTM.crude.local, "VTM.crude.l.local" = VTM.crude.l.local,
                 "VTM.crude.u.local" = VTM.crude.u.local)
  return(output)
}

#-------------------------------------------------------------------------------

# Function to add the constraints at the unit level consequence of both the scenario
# and the fact the row or the column sums zero.

constraints_zeros_local <- function(pjk.crude.local, pjk.crude.l.local, pjk.crude.u.local,
                                       scenario, J0, K0, x, y){
  J <- dim(pjk.crude.local)[1L]
  K <- dim(pjk.crude.local)[2L]
  # zeros constraints due to the scenario
  if((scenario == "raw" & J > J0 & K > K0) |
     (scenario == "regular" & J == J0 & K > K0) |
     (scenario == "ordinary" & J > J0 & K == K0) |
     (scenario == "enriched" & J == J0 & K == K0) |
     (scenario == "semifull")) {
    pjk.crude.local[J, K, ] <- pjk.crude.l.local[J, K, ] <- pjk.crude.u.local[J, K, ] <- 0L
  } else if ((scenario == "regular" & J > J0 & K > K0) |
             (scenario == "enriched" & J > J0 & K == K0) |
             (scenario == "full")) {
    pjk.crude.local[J, K, ] <- pjk.crude.l.local[J, K, ] <- pjk.crude.u.local[J, K, ] <- 0L
    pjk.crude.local[J - 1L, K, ] <- pjk.crude.l.local[J - 1L, K, ] <- pjk.crude.u.local[J - 1L, K, ] <- 0L
  } else if ((scenario == "ordinary" & J > J0 & K > K0) |
             (scenario == "enriched" & J == J0 & K > K0)) {
    pjk.crude.local[J, K, ] <- pjk.crude.l.local[J, K, ] <- pjk.crude.u.local[J, K, ] <- 0L
    pjk.crude.local[J, K - 1L, ] <- pjk.crude.l.local[J, K - 1L, ] <- pjk.crude.u.local[J, K -1L, ] <- 0L
  } else if ((scenario == "enriched" & J > J0 & K > K0) |
             (scenario == "gold")) {
    pjk.crude.local[J, K, ] <- pjk.crude.l.local[J, K, ] <- pjk.crude.u.local[J, K, ] <- 0L
    pjk.crude.local[J, K - 1L, ] <- pjk.crude.l.local[J, K - 1L, ] <- pjk.crude.u.local[J, K -1L, ] <- 0L
    pjk.crude.local[J - 1L, K, ] <- pjk.crude.l.local[J - 1L, K, ] <- pjk.crude.u.local[J - 1L, K, ] <- 0L
    pjk.crude.local[J - 1L, K - 1L, ] <- pjk.crude.l.local[J - 1L, K - 1L, ] <- pjk.crude.u.local[J - 1L, K - 1L, ] <- 0L
  }
  # zeros constraints due to nulls of row or column aggregates
  zeros.row <- which(x == 0L) %% nrow(x)
  zeros.col <- which(y == 0L) %% nrow(y)
  for (ii in zeros.row){
    jj <- which(x[ii, ] == 0L)
    pjk.crude.local[jj, , ii] <- pjk.crude.l.local[jj, , ii] <- pjk.crude.u.local[jj, , ii] <- 0L
  }
  for (ii in zeros.col){
    kk <- which(y[ii, ] == 0L)
    pjk.crude.local[, kk, ii] <- pjk.crude.l.local[, kk, ii] <- pjk.crude.u.local[, kk, ii] <- 0L
  }

  output <- list("pjk.crude.local" = pjk.crude.local, "pjk.crude.l.local" = pjk.crude.l.local,
                 "pjk.crude.u.local" = pjk.crude.u.local)
  return(output)
}

#-------------------------------------------------------------------------------

# Function for extract a sample of a crude transfer matrix
extract_sample <- function(TM.low, TM.upp, B){
  J <- nrow(TM.low); K <- ncol(TM.low)
  output <- array(NA, dim = c(J, K, B))
  for (j in 1L:J){
    for (k in 1L:K){
      output[j, k, ] <- stats::runif(B, min = TM.low[j, k], max = TM.upp[j, k])
    }
  }
  return(output)
}

#-------------------------------------------------------------------------------

# Function for estimating the confidence intervals of the adjusted transitions probabilities
interval_transitions <- function(muestra, vector.fila, vector.columna, tol, confidence){
  B <- dim(muestra)[3]
  for (b in 1L:B){
    TM <- IPF(muestra[, , b], vector.columna, vector.fila, precision = tol)
    muestra[, , b] <- TM/rowSums(TM)
  }
  TM.upp <- apply(muestra, c(1L, 2L), stats::quantile, probs = 1L - (1L - confidence)/2L)
  TM.low <- apply(muestra, c(1L ,2L), stats::quantile, probs = (1L - confidence)/2L)
  output <- list("TM.low" = TM.low, "TM.upp" = TM.upp)
  return(output)
}
#-------------------------------------------------------------------------------

# Function for estimating the confidence intervals of the adjusted transfers votes
interval_transfers <- function(muestra, vector.fila, vector.columna, tol, confidence){
  B <- dim(muestra)[3]
  for (b in 1L:B){
    TM <- IPF(muestra[, , b], vector.columna, vector.fila, precision = tol)
  }
  TM.upp <- apply(TM, c(1L, 2L), stats::quantile, probs = 1L - (1L - confidence)/2L)
  TM.low <- apply(TM, c(1L ,2L), stats::quantile, probs = (1L - confidence)/2L)
  output <- list("TM.low" = TM.low, "TM.upp" = TM.upp)
  return(output)
}
#-------------------------------------------------------------------------------

# Function to estimate the I RxC local transition matrices using logit transformations

Thomsen_local_logit_2x2 <- function(votes1, total1, votes2, total2, confidence = 0.95,
                                    Yule.aprox = FALSE){
  # Inputs:
  # votes1: matrix/vector nx1 with the votes in the first election of the party of interest
  # total1: matrix/vector nx1 with the total votes recorded in the first election
  # votes2: matrix/vector nx1 with the votes in the second election of the party of interest
  # total2: matrix/vector nx1 with the total votes recorded in the second election
  # confidence: level of confidence for the interval estimates. If NULL: not computed.
  # Yule.aprox = TRUE/FALSE argument indicating if either formula (3.44) or (3.46) of Thomsen (1987) should be use to estimate cross-proportions. Default FALSE, formula (3.44).
  # The variables should be in votes (counts)
  # This function estimates the transfers between party1 and party2 in the I units and in the whole space
  # Output: A list with the following elements
  # PTM: A 2x2 matrix with the transitions probabilities between votes and total
  # PTM_local: A 2x2xI array with the transitions probabilities between votes and total in the I units
  # PTM_low: A 2x2 matrix with the upper bounds of transitions probabilities between votes and total
  # PTM_low_local: A 2x2xI matrix with the upper bounds of transitions probabilities between votes and total in the I units
  # PTM_high: A 2x2 matrix with the lower bounds of transitions probabilities between votes and total
  # PTM_high_local: A 2x2xI matrix with the lower bounds of transitions probabilities between votes and total in the I units
  # pjk: A 2x2xI array with the punctual estimates of cross-proportions of simultaneously voting for party 1 and party 2.
  # pjk_low: A 2x2xI vector with the minimum punctual estimates the cross-proportions of simultaneously voting for party 1 and party 2.
  # pjk_high: A 2x2xI vector with the maximum punctual estimates of cross-proportions of simultaneously voting for party 1 and party 2.
  # cor: An estimate of correlation between p(votes1) and p(votes2) in the transformed scale.

  #  This function is partially based on the .m (MATLAB) and .ado (STATA)
  #  functions written by Won-ho Park, University of Michigan in July-Aug. 2002
  votes1 <- as.vector(votes1)
  votes2 <- as.vector(votes2)
  total1 <- as.vector(total1)
  total2 <- as.vector(total1)

  I <- length(votes2)
  weight <- total2/sum(total2)
  t1 <- sum(votes1)/sum(total1)
  t2 <- sum(votes2)/sum(total2)

  # To avoid taking logs from zero values
  votes1[votes1 == 0L] <- 0.5
  votes2[votes2 == 0L] <- 0.5
  total1[total1 == votes1] <- total1[total1 == votes1] + 0.5
  total2[total2 == votes2] <- total2[total2 == votes2] + 0.5

  # Logit transformation of the proportions
  lv1 <- log(votes1/(total1 - votes1))
  lv2 <- log(votes2/(total2 - votes2))
  meanlv1 <- sum(weight * lv1)
  meanlv2 <- sum(weight * lv2)
  varlv1 <- sum(I*weight * (lv1 - meanlv1)^2L)/(I - 1L)
  varlv2 <- sum(I*weight * (lv2 - meanlv2)^2L)/(I - 1L)

  # Pearson correlation of logit transformation
  r <- sum(weight * (lv1 - meanlv1) * (lv2 - meanlv2))/sqrt(varlv1*varlv2)

  # Whole space votes
  # Integration binormal
  core <- pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r)
  if (Yule.aprox){
    # Tetrachoric formula
    k <- sqrt((1L + 2L*r*t1 + 2L*r*t2 - r)^2L - 8L*r*(1L + r)*t1*t2)
    core <- (1L + 2L*r*t1 + 2L*r*t2 - r - k)/(4L * r)
  }
  p2p <-  as.numeric(core)/t1  # Estimate of vote transfer from party 1 to party 2
  t2p <- (t2 - core)/(1L - t1) # Estimate of vote transfer from total(rest) 1 to party 2
  trans.matrix <- matrix(c(p2p, t2p, 1L - p2p, 1L - t2p), 2L, 2L)
  trans.matrix_low <- trans.matrix_high <- array(NA, c(2L, 2L))

  # Units
  trans.matrix_local <- array(NA, c(2L, 2L, I))
  trans.matrix_low_local <- trans.matrix_high_local <- array(NA, c(2L, 2L, I))
  pjk <- pjk_low <- pjk_high <- array(NA, c(2L, 2L, I))

  for (ii in 1L:I){
    t1 <- votes1[ii]/total1[ii]
    t2 <- votes2[ii]/total2[ii]
    # Integration binormal
    core <- pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r)
    if (Yule.aprox){
      # Tetrachoric formula
      k <- sqrt((1L + 2L*r*t1 + 2L*r*t2 - r)^2L - 8L*r*(1L + r)*t1*t2)
      core <- (1L + 2L*r*t1 + 2L*r*t2 - r - k)/(4L * r)
    }
    # p12[ii] <- as.numeric(core)  # Estimate of cross-probability from party 1 to party 2
    p2p <-  as.numeric(core)/t1  # Estimate of proportion transfer from party 1 to party 2
    t2p <- (t2 - core)/(1L - t1) # Estimate of proportion transfer from total(rest) 1 to party 2
    trans.matrix_local[, , ii] <- matrix(c(p2p, t2p, 1L - p2p, 1L - t2p), 2, 2)
  }

  # Standard deviations
  if (!is.null(confidence)) {
    z <- stats::qnorm(1L - (1L - confidence)/2L)
    # z_hi <- log((1L + r)/(1L - r))/2L + z/sqrt(I - 3L)
    # z_low <- log((1L + r)/(1 - r))/2L - z/sqrt(I - 3L)
    # r_hi <- (exp(2L * z_hi) - 1L)/(exp(2L * z_hi) + 1L)
    # r_low <- (exp(2L * z_low) - 1L)/(exp(2L * z_low) +1L)
    r_hi <- tanh(log((1L + r)/(1L - r))/2L + z/sqrt(I - 3L))
    r_low <- tanh(log((1L + r)/(1L - r))/2L - z/sqrt(I - 3L))

    # Integration binormal
    p2p_hi <- as.numeric(pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r_hi))/t1
    p2p_low <- as.numeric(pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r_low))/t1
    if (Yule.aprox){
      # Tetrachoric formula
      k <- sqrt((1L + 2L*r_hi*t1 + 2L*r_hi*t2 - r_hi)^2L - 8L*r_hi*(1L + r_hi)*t1*t2)
      p2p_hi <- (1L + 2L*r_hi*t1 + 2L*r_hi*t2 - r_hi - k)/(4L * r_hi)/t1
      k <- sqrt((1L + 2L*r_low*t1 + 2L*r_low*t2 - r_low)^2L - 8L*r_low*(1L + r_low)*t1*t2)
      p2p_low <- (1L + 2L*r_low*t1 + 2L*r_low*t2 - r_low - k)/(4L * r_low)/t1
    }
    t2p_hi <- (t2 - p2p_low*t1)/(1L - t1)
    t2p_low <- (t2 - p2p_hi*t1)/(1L - t1)

    trans.matrix_low <- matrix(c(p2p_low, t2p_low, 1L - p2p_hi, 1L - t2p_hi), 2L, 2L)
    trans.matrix_high <- matrix(c(p2p_hi, t2p_hi, 1L - p2p_low, 1L - t2p_low), 2L, 2L)
    # By units

    for (ii in 1L:I){
      t1 <- votes1[ii]/total1[ii]
      t2 <- votes2[ii]/total2[ii]

      # Integration binormal
      p2p_hi <- as.numeric(pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r_hi))/t1
      p2p_low <- as.numeric(pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r_low))/t1
      if (Yule.aprox){
        # Tetrachoric formula
        k <- sqrt((1L + 2L*r_hi*t1 + 2L*r_hi*t2 - r_hi)^2L - 8L*r_hi*(1L + r_hi)*t1*t2)
        p2p_hi <- (1L + 2L*r_hi*t1 + 2L*r_hi*t2 - r_hi - k)/(4L * r_hi)/t1
        k <- sqrt((1L + 2L*r_low*t1 + 2L*r_low*t2 - r_low)^2L - 8L*r_low*(1L + r_low)*t1*t2)
        p2p_low <- (1L + 2L*r_low*t1 + 2L*r_low*t2 - r_low - k)/(4L * r_low)/t1
      }

      t2p_hi <- (t2 - p2p_low*t1)/(1L - t1)
      t2p_low <- (t2 - p2p_hi*t1)/(1L - t1)

      trans.matrix_low_local[, , ii] <- matrix(c(p2p_low, t2p_low, 1L - p2p_hi, 1L - t2p_hi), 2L, 2L)
      trans.matrix_high_local[, , ii] <- matrix(c(p2p_hi, t2p_hi, 1L - p2p_low, 1L - t2p_low), 2L, 2L)
    }

  }

  total.rows <- array(as.vector(kronecker(c(1L ,1L), rbind(votes1/total1, 1L - votes1/total1))),
                      c(2L, 2L, I))
  pjk <- trans.matrix_local*total.rows
  pjk_low <- trans.matrix_low_local*total.rows
  pjk_high <- trans.matrix_high_local*total.rows

  output <- list("PTM" = trans.matrix, "PTM_low" = trans.matrix_low,
                 "PTM_high" = trans.matrix_high, "PTM_local" = trans.matrix_local,
                 "PTM_low_local" = trans.matrix_low_local, "PTM_high_local" = trans.matrix_high_local,
                 "pjk" = pjk, "pjk_low" = pjk_low, "pjk_high" = pjk_high, "cor" = r)

  return(output)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to estimate the I RxC local transition matrices using probit transformations

Thomsen_local_probit_2x2 <- function(votes1, total1, votes2, total2, confidence = 0.95,
                                     Yule.aprox = FALSE){
  # Inputs:
  # votes1: matrix/vector nx1 with the votes in the first election of the party of interest
  # total1: matrix/vector nx1 with the total votes recorded in the first election
  # votes2: matrix/vector nx1 with the votes in the second election of the party of interest
  # total2: matrix/vector nx1 with the total votes recorded in the second election
  # confidence: level of confidence for the interval estimates. If NULL: not computed.
  # Yule.aprox = TRUE/FALSE argument indicating if either formula (3.44) or (3.46) of Thomsen (1987) should be use to estimate cross-proportions. Default FALSE, formula (3.44).
  # The variables should be in votes (counts)
  # This function estimates the transfers between party1 and party2 in the I units and in the whole space
  # Output: A list with the following elements
  # PTM: A 2x2 matrix with the transitions probabilities between votes and total
  # PTM_local: A 2x2xI array with the transitions probabilities between votes and total in the I units
  # PTM_low: A 2x2 matrix with the upper bounds of transitions probabilities between votes and total
  # PTM_low_local: A 2x2xI matrix with the upper bounds of transitions probabilities between votes and total in the I units
  # PTM_high: A 2x2 matrix with the lower bounds of transitions probabilities between votes and total
  # PTM_high_local: A 2x2xI matrix with the lower bounds of transitions probabilities between votes and total in the I units
  # pjk: A 2x2xI array with the punctual estimates of cross-proportions of simultaneously voting for party 1 and party 2.
  # pjk_low: A 2x2xI vector with the minimum punctual estimates the cross-proportions of simultaneously voting for party 1 and party 2.
  # pjk_high: A 2x2xI vector with the maximum punctual estimates of cross-proportions of simultaneously voting for party 1 and party 2.
  # cor: An estimate of correlation between p(votes1) and p(votes2) in the transformed scale.

  #  This function is partially based on the .m (MATLAB) and .ado (STATA)
  #  functions written by Won-ho Park, University of Michigan in July-Aug. 2002
  votes1 <- as.vector(votes1)
  votes2 <- as.vector(votes2)
  total1 <- as.vector(total1)
  total2 <- as.vector(total1)

  I <- length(votes2)
  weight <- total2/sum(total2)
  t1 <- sum(votes1)/sum(total1)
  t2 <- sum(votes2)/sum(total2)

  # To avoid taking probits from one or zero values
  votes1[votes1 == 0L] <- 0.5
  votes2[votes2 == 0L] <- 0.5
  total1[total1 == votes1] <- total1[total1 == votes1] + 0.5
  total2[total2 == votes2] <- total2[total2 == votes2] + 0.5

  # Probit transformation of the proportions
  pv1 <- stats::qnorm(votes1/total1)
  pv2 <- stats::qnorm(votes2/total2)
  meanpv1 <- sum(weight * pv1)
  meanpv2 <- sum(weight * pv2)
  varpv1 <- sum(I*weight * (pv1 - meanpv1)^2L)/(I - 1L)
  varpv2 <- sum(I*weight * (pv2 - meanpv2)^2L)/(I - 1L)

  # Pearson correlation of probit transformations
  r <- sum(weight * (pv1 - meanpv1) * (pv2 - meanpv2))/sqrt(varpv1*varpv2)

  # Whole space
  # Integration binormal
  core <- pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r)
  if (Yule.aprox){
    # Tetrachoric formula
    k <- sqrt((1L + 2L*r*t1 + 2L*r*t2 - r)^2L - 8L*r*(1L + r)*t1*t2)
    core <- (1L + 2L*r*t1 + 2L*r*t2 - r - k)/(4L * r)
  }
  p2p <-  as.numeric(core)/t1  # Estimate of proportion transfer from party 1 to party 2
  t2p <- (t2 - core)/(1L - t1) # Estimate of proportion transfer from total(rest) 1 to party 2
  trans.matrix <- matrix(c(p2p, t2p, 1L - p2p, 1L - t2p), 2L, 2L)
  trans.matrix_low <- trans.matrix_high <- array(NA, c(2L, 2L))

  # Units
  trans.matrix_local <- array(NA, c(2L, 2L, I))
  trans.matrix_low_local <- trans.matrix_high_local <- array(NA, c(2L, 2L, I))
  pjk <- pjk_low <- pjk_high <- array(NA, c(2L, 2L, I))

  for (ii in 1L:I){
    t1 <- votes1[ii]/total1[ii]
    t2 <- votes2[ii]/total2[ii]
    core <- pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r)
    if (Yule.aprox){
      # Tetrachoric formula
      k <- sqrt((1L + 2L*r*t1 + 2L*r*t2 - r)^2L - 8L*r*(1L + r)*t1*t2)
      core <- (1L + 2L*r*t1 + 2L*r*t2 - r - k)/(4L * r)
    }
    # p12[ii] <- as.numeric(core)  # Estimate of cross-probability from party 1 to party 2
    p2p <-  as.numeric(core)/t1  # Estimate of proportion transfer from party 1 to party 2
    t2p <- (t2 - core)/(1L - t1) # Estimate of proportion transfer from total(rest) 1 to party 2
    trans.matrix_local[, , ii] <- matrix(c(p2p, t2p, 1L - p2p, 1L - t2p), 2, 2)
  }

  # Standard deviations
  if (!is.null(confidence)) {

    z <- stats::qnorm(1L - (1L - confidence)/2L)
    r_hi <- tanh(log((1L + r)/(1L - r))/2L + z/sqrt(I - 3L))
    r_low <- tanh(log((1L + r)/(1L - r))/2L - z/sqrt(I - 3L))

    # Whole space
    # Integration binormal
    p2p_hi <- as.numeric(pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r_hi))/t1
    p2p_low <- as.numeric(pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r_low))/t1
    if (Yule.aprox){
      # Tetrachoric formula
      k <- sqrt((1L + 2L*r_hi*t1 + 2L*r_hi*t2 - r_hi)^2L - 8L*r_hi*(1L + r_hi)*t1*t2)
      p2p_hi <- (1L + 2L*r_hi*t1 + 2L*r_hi*t2 - r_hi - k)/(4L * r_hi)/t1
      k <- sqrt((1L + 2L*r_low*t1 + 2L*r_low*t2 - r_low)^2L - 8L*r_low*(1L + r_low)*t1*t2)
      p2p_low <- (1L + 2L*r_low*t1 + 2L*r_low*t2 - r_low - k)/(4L * r_low)/t1
    }
    t2p_hi <- (t2 - p2p_low*t1)/(1L - t1)
    t2p_low <- (t2 - p2p_hi*t1)/(1L - t1)

    trans.matrix_low <- matrix(c(p2p_low, t2p_low, 1L - p2p_hi, 1L - t2p_hi), 2L, 2L)
    trans.matrix_high <- matrix(c(p2p_hi, t2p_hi, 1L - p2p_low, 1L - t2p_low), 2L, 2L)

    # By units
    trans.matrix_low_local <- trans.matrix_high_local <- array(NA, c(2L, 2L, I))
    for (ii in 1L:I){
      t1 <- votes1[ii]/total1[ii]
      t2 <- votes2[ii]/total2[ii]
      # Integration binormal
      p2p_hi <- as.numeric(pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r_hi))/t1
      p2p_low <- as.numeric(pnorm2d(stats::qnorm(t1), stats::qnorm(t2), rho = r_low))/t1
      if (Yule.aprox){
        # Tetrachoric formula
        k <- sqrt((1L + 2L*r_hi*t1 + 2L*r_hi*t2 - r_hi)^2L - 8L*r_hi*(1L + r_hi)*t1*t2)
        p2p_hi <- (1L + 2L*r_hi*t1 + 2L*r_hi*t2 - r_hi - k)/(4L * r_hi)/t1
        k <- sqrt((1L + 2L*r_low*t1 + 2L*r_low*t2 - r_low)^2L - 8L*r_low*(1L + r_low)*t1*t2)
        p2p_low <- (1L + 2L*r_low*t1 + 2L*r_low*t2 - r_low - k)/(4L * r_low)/t1
      }

      t2p_hi <- (t2 - p2p_low*t1)/(1L - t1)
      t2p_low <- (t2 - p2p_hi*t1)/(1L - t1)

      trans.matrix_low_local[, , ii] <- matrix(c(p2p_low, t2p_low, 1L - p2p_hi, 1L - t2p_hi), 2L, 2L)
      trans.matrix_high_local[, , ii] <- matrix(c(p2p_hi, t2p_hi, 1L - p2p_low, 1L - t2p_low), 2L, 2L)
    }

  }

  total.rows <- array(as.vector(kronecker(c(1L ,1L), rbind(votes1/total1, 1L - votes1/total1))),
                      c(2L, 2L, I))
  pjk <- trans.matrix_local*total.rows
  pjk_low <- trans.matrix_low_local*total.rows
  pjk_high <- trans.matrix_high_local*total.rows

  output <- list("PTM" = trans.matrix, "PTM_low" = trans.matrix_low,
                 "PTM_high" = trans.matrix_high, "PTM_local" = trans.matrix_local,
                 "PTM_low_local" = trans.matrix_low_local, "PTM_high_local" = trans.matrix_high_local,
                 "pjk" = pjk, "pjk_low" = pjk_low, "pjk_high" = pjk_high, "cor" = r)

  return(output)
}


#-------------------------------------------------------------------------------
# Function to reweighing the cross-proportions pjk in the I units given a row and a column
# party as references using the approach proposed by Thomsen (1987) through equations (4.17)-(4.22)
# but without approximation, allowing both transformations: logit and probit.

reweighting_pjk <- function(pjk.array, ref1, ref2, scale, x, y, tol = 10^-6){
  I <- dim(pjk.array)[3L]
  J <- dim(pjk.array)[1L]
  K <- dim(pjk.array)[2L]

  parties1 <- c(1L:J)[-ref1]
  parties2 <- c(1L:K)[-ref2]

  p.obs.rows <- x/rowSums(x)
  p.obs.cols <- y/rowSums(y)
  weight <- rowSums(x)/sum(x)

  suma0 <- 0L
  suma1 <- sum(pjk.array)
  dif <- abs(suma0 - suma1)
  cont <- 0L

  while(dif > tol & cont < 10000L){
    if (scale == "logit"){
      for (jj in parties1){
        for (kk in parties2){
          pjk.temp <- pjk.array[c(jj, ref1), c(kk, ref2), ]
          suma.rows <- apply(pjk.temp, c(1L, 3L), sum)
          suma.cols <- apply(pjk.temp, c(2L, 3L), sum)
          ps1 <- suma.rows[1L, ]/suma.rows[2L, ]
          ps2 <- suma.cols[1L, ]/suma.cols[2L, ]
          weight <- rowSums(x) * apply(suma.rows, 2, sum)
          weight <- weight/sum(weight)
          valid <- abs(ps1) != 0L & !is.infinite(ps1) & ps2 != 0L & !is.infinite(ps2) & !is.na(ps1) & !is.na(ps2)
          valid <- valid & suma.rows[1L, ] >= 1L & suma.rows[1L, ] == 0L & suma.cols[1L, ] >= 1L & suma.cols[1L, ] == 0L
          logit.rows <- log(ps1[valid])
          logit.cols <- log(ps2[valid])
          mean1 <- sum(weight[valid] * logit.rows)/sum(weight[valid])
          mean2 <- sum(weight[valid] * logit.cols)/sum(weight[valid])
          var1 <- sum(weight[valid] * (logit.rows - mean1)^2L)/sum(weight[valid])
          var2 <- sum(weight[valid] * (logit.cols - mean2)^2L)/sum(weight[valid])
          # Pearson correlation of logit transformations
          cor.e <- sum(weight[valid] * (logit.rows - mean1) * (logit.cols - mean2))/sqrt(var1*var2)
          no.0 <- valid & pjk.array[jj, kk, ] != 0L
          pjk.array[jj, kk, no.0] <- pnorm2d(stats::qnorm(suma.rows[1L, no.0]),
                                             stats::qnorm(suma.cols[1L, no.0]), rho = cor.e)
        }
      }
    } else {
      for (jj in parties1){
        for (kk in parties2){
          pjk.temp <- pjk.array[c(jj, ref1), c(kk, ref2), ]
          suma.rows <- apply(pjk.temp, c(1L, 3L), sum)
          suma.cols <- apply(pjk.temp, c(2L, 3L), sum)
          ps1 <- suma.rows[1L, ]/apply(suma.rows, 2, sum)
          ps2 <- suma.cols[1L, ]/apply(suma.cols, 2, sum)
          weight <- rowSums(x) * apply(suma.rows, 2, sum)
          weight <- weight/sum(weight)
          valid <- ps1 != 0L & !is.infinite(ps1) & ps2 != 0L & !is.infinite(ps2) & ps1 != 1L & ps2 != 1L & !is.na(ps1) & !is.na(ps2)
          valid <- valid & suma.rows[1L, ] >= 1L & suma.rows[1L, ] == 0L & suma.cols[1L, ] >= 1L & suma.cols[1L, ] == 0L
          probit.rows <- stats::qnorm(ps1[valid])
          probit.cols <- stats::qnorm(ps2[valid])
          mean1 <- sum(weight[valid] * probit.rows)/sum(weight[valid])
          mean2 <- sum(weight[valid] * probit.cols)/sum(weight[valid])
          #var1 <- sum(I*weight[valid] * (probit.rows - mean1)^2L)/(I - 1L)
          #var2 <- sum(I*weight[valid] * (probit.cols - mean2)^2L)/(I - 1L)
          var1 <- sum(weight[valid] * (probit.rows - mean1)^2L)/sum(weight[valid])
          var2 <- sum(weight[valid] * (probit.cols - mean2)^2L)/sum(weight[valid])
          # Pearson correlation of probit transformations
          cor.e <- sum(weight[valid] * (probit.rows - mean1) * (probit.cols - mean2))/sqrt(var1*var2)
          no.0 <- valid & pjk.array[jj, kk, ] != 0L
          pjk.array[jj, kk, no.0] <- pnorm2d(stats::qnorm(suma.rows[1L, no.0]),
                                             stats::qnorm(suma.cols[1L, no.0]), rho = cor.e)
        }
      }
    }

    pjk.array[pjk.array < 0] <- 0L
    p.est.rows <- t(apply(pjk.array, c(1L, 3L), sum))
    p.est.cols <- t(apply(pjk.array, c(2L, 3L), sum))

    for (kk in parties2){
      no.0 <- p.est.cols[, kk] != 0L
      pjk.array[ref1, kk, no.0] <- pjk.array[ref1, kk, no.0]*p.obs.cols[no.0, kk]/p.est.cols[no.0, kk]
    }
    for (jj in parties1){
      no.0 <- p.est.rows[, jj] != 0L
      pjk.array[jj, ref2, no.0] <- pjk.array[jj, ref2, no.0]*p.obs.rows[no.0, jj]/p.est.rows[no.0, jj]
    }
    cont <- cont + 1L
    suma0 <- suma1
    suma1 <- sum(pjk.array)
    dif <- abs(suma0 - suma1)
  }

  for(ii in 1L:I){
     pjk.array[, , ii] <- IPF(pjk.array[, , ii]*sum(y[ii, ]), y[ii, ], x[ii, ], tol)
     dif.r <- sum(abs(x[ii, ] - rowSums(pjk.array[, , ii])))
     dif.c <- sum(abs(y[ii, ] - colSums(pjk.array[, , ii])))
     if (max(dif.r, dif.c) > 0.01){
       pjk.array[, , ii] <- IPF2(pjk.array[, , ii], y[ii, ], x[ii, ])
     }
  }

  output <- list("vjk.array" = pjk.array, "iter" = cont)
  return(output)

}


#-------------------------------------------------------------------------------
# Function to reweighing the cross-proportions pjk in the I units given a row and a column
# party as references using the approach proposed by Thomsen (1987) through equations (4.17)-(4.22)
# allowing both transformations: logit and probit.

reweighting_pjk_Yule <- function(pjk.array, ref1, ref2, scale, x, y, tol = 10^-6){
  I <- dim(pjk.array)[3L]
  J <- dim(pjk.array)[1L]
  K <- dim(pjk.array)[2L]

  parties1 <- c(1L:J)[-ref1]
  parties2 <- c(1L:K)[-ref2]

  p.obs.rows <- x/rowSums(x)
  p.obs.cols <- y/rowSums(y)
  weight <- rowSums(x)/sum(x)

  suma0 <- 0L
  suma1 <- sum(pjk.array)
  dif <- abs(suma0 - suma1)
  cont <- 0L

  while(dif > tol & cont < 10000L){
    if (scale == "logit"){
      for (jj in parties1){
        for (kk in parties2){
          pjk.temp <- pjk.array[c(jj, ref1), c(kk, ref2), ]
          suma.rows <- apply(pjk.temp, c(1L, 3L), sum)
          suma.cols <- apply(pjk.temp, c(2L, 3L), sum)
          ps1 <- suma.rows[1L, ]/suma.rows[2L, ]
          ps2 <- suma.cols[1L, ]/suma.cols[2L, ]
          weight <- rowSums(x) * apply(suma.rows, 2, sum)
          weight <- weight/sum(weight)
          valid <- abs(ps1) != 0L & !is.infinite(ps1) & ps2 != 0L & !is.infinite(ps2) & !is.na(ps1) & !is.na(ps2)
          valid <- valid & suma.rows[1L, ] >= 1L & suma.rows[1L, ] == 0L & suma.cols[1L, ] >= 1L & suma.cols[1L, ] == 0L
          logit.rows <- log(ps1[valid])
          logit.cols <- log(ps2[valid])
          mean1 <- sum(weight[valid] * logit.rows)/sum(weight[valid])
          mean2 <- sum(weight[valid] * logit.cols)/sum(weight[valid])
          var1 <- sum(weight[valid] * (logit.rows - mean1)^2L)/sum(weight[valid])
          var2 <- sum(weight[valid] * (logit.cols - mean2)^2L)/sum(weight[valid])
          # Pearson correlation of logit transformations
          cor.e <- sum(weight[valid] * (logit.rows - mean1) * (logit.cols - mean2))/sqrt(var1*var2)
          no.0 <- valid & pjk.array[jj, kk, ] != 0L
          k <- (1L + cor.e)/(1L - cor.e)
          pjk.array[jj, kk, no.0] <- k*(pjk.array[jj, ref2, no.0]*
                                        pjk.array[ref1, kk, no.0])/pjk.array[ref1, ref2, no.0]
        }
      }
    } else {
      for (jj in parties1){
        for (kk in parties2){
          pjk.temp <- pjk.array[c(jj, ref1), c(kk, ref2), ]
          suma.rows <- apply(pjk.temp, c(1L, 3L), sum)
          suma.cols <- apply(pjk.temp, c(2L, 3L), sum)
          ps1 <- suma.rows[1L, ]/apply(suma.rows, 2, sum)
          ps2 <- suma.cols[1L, ]/apply(suma.cols, 2, sum)
          weight <- rowSums(x) * apply(suma.rows, 2, sum)
          weight <- weight/sum(weight)
          valid <- ps1 != 0L & !is.infinite(ps1) & ps2 != 0L & !is.infinite(ps2) & ps1 != 1L & ps2 != 1L & !is.na(ps1) & !is.na(ps2)
          valid <- valid & suma.rows[1L, ] >= 1L & suma.rows[1L, ] == 0L & suma.cols[1L, ] >= 1L & suma.cols[1L, ] == 0L
          probit.rows <- stats::qnorm(ps1[valid])
          probit.cols <- stats::qnorm(ps2[valid])
          mean1 <- sum(weight[valid] * probit.rows)/sum(weight[valid])
          mean2 <- sum(weight[valid] * probit.cols)/sum(weight[valid])
          #var1 <- sum(I*weight[valid] * (probit.rows - mean1)^2L)/(I - 1L)
          #var2 <- sum(I*weight[valid] * (probit.cols - mean2)^2L)/(I - 1L)
          var1 <- sum(weight[valid] * (probit.rows - mean1)^2L)/sum(weight[valid])
          var2 <- sum(weight[valid] * (probit.cols - mean2)^2L)/sum(weight[valid])
          # Pearson correlation of probit transformations
          cor.e <- sum(weight[valid] * (probit.rows - mean1) * (probit.cols - mean2))/sqrt(var1*var2)
          no.0 <- valid & pjk.array[jj, kk, ] != 0L
          k <- (1L + cor.e)/(1L - cor.e)
          pjk.array[jj, kk, no.0] <- k*(pjk.array[jj, ref2, no.0]*
                                        pjk.array[ref1, kk, no.0])/pjk.array[ref1, ref2, no.0]
        }
      }
    }

    pjk.array[pjk.array < 0] <- 0L
    p.est.rows <- t(apply(pjk.array, c(1L, 3L), sum))
    p.est.cols <- t(apply(pjk.array, c(2L, 3L), sum))

    for (kk in parties2){
      no.0 <- p.est.cols[, kk] != 0L
      pjk.array[ref1, kk, no.0] <- pjk.array[ref1, kk, no.0]*p.obs.cols[no.0, kk]/p.est.cols[no.0, kk]
    }
    for (jj in parties1){
      no.0 <- p.est.rows[, jj] != 0L
      pjk.array[jj, ref2, no.0] <- pjk.array[jj, ref2, no.0]*p.obs.rows[no.0, jj]/p.est.rows[no.0, jj]
    }
    cont <- cont + 1L
    suma0 <- suma1
    suma1 <- sum(pjk.array)
    dif <- abs(suma0 - suma1)
  }

  for(ii in 1L:I){
    pjk.array[, , ii] <- IPF(pjk.array[, , ii]*sum(y[ii, ]), y[ii, ], x[ii, ], tol)
    dif.r <- sum(abs(x[ii, ] - rowSums(pjk.array[, , ii])))
    dif.c <- sum(abs(y[ii, ] - colSums(pjk.array[, , ii])))
    if (max(dif.r, dif.c) > 0.01){
      pjk.array[, , ii] <- IPF2(pjk.array[, , ii], y[ii, ], x[ii, ])
    }
  }

  output <- list("vjk.array" = pjk.array, "iter" = cont)
  return(output)

}

#-------------------------------------------------------------------------------
# Function to compute the global solution corresponding to each pair of references.

solutions_by_reference <- function(vjk.units.multi, x0, y0){
  J0 <- ncol(x0)
  K0 <- ncol(y0)

  vjk_r1_r2 <- array(NA, dim(vjk.units.multi)[-3])
  for (hh in 1L:(J0*K0)){
    temp <- vjk.units.multi[ , , , hh]
    vjk_r1_r2[, , hh] <- apply(temp, c(1L, 2L), sum)
  }

  return("vjk_r1_r2" = vjk_r1_r2)
}

#-------------------------------------------------------------------------------
# Function to combine/average local/global solutions to reach a global solution using global weights.
average_global_weights <- function(vjk.units.multi, VTM.crude, x0, y0, correlations){
  J0 <- ncol(x0)
  K0 <- ncol(y0)
  J <- nrow(VTM.crude)
  K <- ncol(VTM.crude)

  # Reference Cell Number of Voters: RCNV
  vjk_r1_r2 <- solutions_by_reference(vjk.units.multi, x0, y0)
  weights <- as.vector(t(VTM.crude[1L:J0, 1L:K0] * colSums(x0)))
  weights <- weights/sum(weights)
  W <- array(rep(weights, each = J*K), dim(vjk_r1_r2))
  RCNV <- apply(vjk_r1_r2 * W, c(1L, 2L), sum)

  # Square Root Reference Cell Number of Voters: SQRCNV
  weights <- as.vector(t(sqrt(VTM.crude[1L:J0, 1L:K0] * colSums(x0))))
  weights <- weights/sum(weights)
  W <- array(rep(weights, each = J*K), dim(vjk_r1_r2))
  SQRCNV <- apply(vjk_r1_r2 * W, c(1L, 2L), sum)

  # Square Root Reference Margins: SQRM
  weights <- NULL
  for (jj in 1L:J0){
    for (kk in 1L:K0)
    weights <- c(weights, sqrt(colSums(x0)[jj] * colSums(y0)[kk]))
  }
  weights <- weights/sum(weights)
  W <- array(rep(weights, each = J*K), dim(vjk_r1_r2))
  SQRM <- apply(vjk_r1_r2 * W, c(1L, 2L), sum)

  # Correlation Reference: CR
  weights <- abs(as.vector(t(correlations[1L:J0, 1L:K0])))
  weights <- weights/sum(weights)
  W <- array(rep(weights, each = J*K), dim(vjk_r1_r2))
  CR <- apply(vjk_r1_r2 * W, c(1L, 2L), sum)

  output <- list("RCNV" = RCNV, "SQRCNV" = SQRCNV, "SQRM" = SQRM, "AVCR" = CR)
  return(output)
}

#-------------------------------------------------------------------------------
# Function to combine/average local solutions to reach a global solution using local weights.
average_local_weights <- function(vjk.units.multi, pjk.crude.local, x0, y0){
  J0 <- ncol(x0)
  K0 <- ncol(y0)
  J <- dim(pjk.crude.local)[1L]
  K <- dim(pjk.crude.local)[2L]

  # Local Reference Cell Number of Voters: LRCNV
  vjk.crude.local <- pjk.crude.local[1L:J0, 1L:K0, ] * array(rep(rowSums(x0), each = J0*K0),
                                                             dim(pjk.crude.local[1L:J0, 1L:K0, ]))
  weights <- array(rep(as.vector(aperm(vjk.crude.local, c(3, 2, 1))), each = J*K), dim(vjk.units.multi))
  LRCNV <- apply(vjk.units.multi * weights, c(1L, 2L, 3L), sum)
  weights <- array(rep(as.vector(t(apply(vjk.crude.local, 3L, sum))), each =J*K), dim(LRCNV))
  LRCNV <- apply(LRCNV / weights, c(1L, 2L), sum)

  # Local Reference Cell Number of Voters: LSQRCNV
  weights <- array(sqrt(rep(as.vector(aperm(vjk.crude.local, c(3, 2, 1))), each = J*K)),
                   dim(vjk.units.multi))
  LSQRCNV <- apply(vjk.units.multi * weights, c(1L, 2L, 3L), sum)
  weights <- apply(weights, c(1L, 2L, 3L), sum)
  LSQRCNV <- apply(LSQRCNV / weights, c(1L, 2L), sum)

  # Local Square Root Reference Margins: LSQRM
  weights <- NULL
  for (jj in 1L:J0){
    for (kk in 1L:K0)
      weights <- c(weights, rep(sqrt(x0[, jj] * y0[, kk]), each = J*K))
  }
  weights <- array(weights, dim(vjk.units.multi))
  LSQRM <- apply(vjk.units.multi * weights, c(1L, 2L, 3L), sum)
  weights <- apply(weights, c(1L, 2L, 3L), sum)
  LSQRM <- apply(LSQRM / weights, c(1L, 2L), sum)

  output <- return(list("LRCNV" = LRCNV, "LSQRCNV" = LSQRCNV, "LSQRM" = LSQRM))
  return(output)
}

#-------------------------------------------------------------------------------
# Function to generate the reference.outputs.
reference_outputs <- function(vjk.units.multi, VTM.crude, pjk.crude.local, x0, y0,
                              correlations){
  # vjk.by.reference
  names1 <- rownames(VTM.crude)
  names2 <- colnames(VTM.crude)
  names3 <- rownames(x0)
  J0 <- ncol(x0)
  K0 <- ncol(y0)
  names4 <- paste0(rep(paste0("r", 1L:J0, sep =""), each = K0), "_c", 1L:K0, sep = "")
  vjk.by.reference <- solutions_by_reference(vjk.units.multi, x0, y0)
  dimnames(vjk.by.reference) <- c(list(names1), list(names2), list(names4))

  # vjk.units.by.reference
  dimnames(vjk.units.multi) <- c(list(names1), list(names2), list(names3), list(names4))

  # vjk.averages
  vjk.averages <- array(NA, c(dim(VTM.crude), 8L))
  vjk.averages[, , 1L] <- apply(apply(vjk.units.multi, c(1L, 2L, 3L), mean), c(1L, 2L), sum)
  averages <- average_global_weights(vjk.units.multi = vjk.units.multi, VTM.crude = VTM.crude, x0 = x0,
                                     y0 = y0, correlations = correlations)
  vjk.averages[, , 2L] <- averages$RCNV
  vjk.averages[, , 3L] <- averages$SQRCNV
  vjk.averages[, , 4L] <- averages$SQRM
  vjk.averages[, , 5L] <- averages$AVCR
  averages <- average_local_weights(vjk.units.multi = vjk.units.multi, pjk.crude.local = pjk.crude.local,
                                    x0 = x0, y0 = y0)
  vjk.averages[, , 6L] <- averages$LRCNV
  vjk.averages[, , 7L] <- averages$LSQRCNV
  vjk.averages[, , 8L] <- averages$LSQRM
  names3 <- c("Mean", "RCNV", "SQRCNV", "SQRM", "AVCR", "LRCNV", "LSQRCNV", "LSQRM")
  dimnames(vjk.averages) <- c(list(names1), list(names2), list(names3))

  return(list("vjk.averages" = vjk.averages, "vjk.by.reference" = vjk.by.reference,
              "vjk.units.by.reference" = vjk.units.multi))
}

#-------------------------------------------------------------------------------
Thomsen_iter_algorithm <- function(pjk.crude.local, Yule.aprox, reference,
                                   scale, x, y, J0, K0, tol){
#  J <- dim(pjk.crude.local)[1]
#  K <- dim(pjk.crude.local)[2]
#  I <- dim(pjk.crude.local)[3]

  if (Yule.aprox){
    reweighting <- reweighting_pjk_Yule
  } else {
    reweighting <- reweighting_pjk
  }

  vjk.units.multi <- array(NA, c(dim(pjk.crude.local), J0*K0))
  iter <- NULL

  if (!is.null(reference)){
    vjk.units <- reweighting(pjk.crude.local, ref1 = reference[1L] , ref2 = reference[2L],
                             scale = scale, x = x, y = y, tol = tol)
    iter <- vjk.units$iter
    vjk.units <- vjk.units$vjk.array
  } else {
    for (jj in 1:J0){
      for(kk in 1:K0){
        vjk.temp <- reweighting(pjk.crude.local, ref1 = jj , ref2 = kk, scale = scale,
                                x = x, y = y, tol = tol)
        vjk.units.multi[, , ,(jj - 1L)*K0 + kk] <- vjk.temp$vjk.array
        iter <- c(iter, vjk.temp$iter)
      }
    }
    vjk.units <- apply(vjk.units.multi, c(1L, 2L, 3L), mean)
 }
  return(list("vjk.units" = vjk.units, "vjk.units.multi" = vjk.units.multi, "iter" = iter))
}

#-------------------------------------------------------------------------------
# Function for extracting samples of (global and local) vjk matrices after applying the Thomsen algorithm
# and combining reference solutions using correlations
extract_muestra_vjk <- function(pjk.low, pjk.upp, Yule.aprox, reference, scale, x, y, J0, K0, tol,
                                correlations = correlations){
  JKI <- prod(dim(pjk.low))
  pjk.crude <- array(stats::runif(JKI, min = as.vector(pjk.low), max = as.vector(pjk.upp)),
                     dim = dim(pjk.low))

  if(is.null(reference)){
    muestra <- Thomsen_iter_algorithm(pjk.crude.local = pjk.crude,
                                      Yule.aprox = Yule.aprox, reference = reference,
                                      scale = scale, x = x, y = y,
                                      J0 = J0, K0 = K0, tol = tol)$vjk.units.multi
    weights <- abs(as.vector(t(correlations[1L:J0, 1L:K0])))
    weights <- weights/sum(weights)
    W <- array(rep(weights, each = JKI), dim(muestra))
    muestra <- apply(muestra * W, c(1L, 2L, 3L), sum)
  } else {
    muestra <- Thomsen_iter_algorithm(pjk.crude.local = pjk.crude,
                                      Yule.aprox = Yule.aprox, reference = reference,
                                      scale = scale, x = x, y = y,
                                      J0 = J0, K0 = K0, tol = tol)$vjk.units
  }
  return("muestra" = muestra)
}

#-------------------------------------------------------------------------------
# Function for creating intervals using the Thomsen approach
# we adopt a conservative approach

intervals_Thomsen <- function(muestra_vjk_local, vjk.units, x, confidence){
# local
  J <- dim(vjk.units)[1]
  K <- dim(vjk.units)[2]
  total.rows <- array(as.vector(kronecker(rep(1L, K), t(x))), dim(vjk.units))
  VTM.l.local <- apply(muestra_vjk_local, c(1L, 2L, 3L),
                       stats::quantile, probs = (1L - confidence)/2L)
  VTM.u.local <- apply(muestra_vjk_local, c(1L, 2L, 3L),
                       stats::quantile, probs = 1L - (1L - confidence)/2L)
  VTM.m.local <- apply(muestra_vjk_local, c(1L, 2L, 3L), mean)
  inf <- VTM.m.local - VTM.l.local
  sup <- VTM.u.local - VTM.m.local
  VTM.l.local <- pmax(pmin(vjk.units - inf, VTM.l.local), 0L)
  VTM.u.local <- pmin(pmax(vjk.units + sup, VTM.u.local), total.rows)
  VTM.l.local <- VTM.l.local/total.rows
  VTM.u.local <- VTM.u.local/total.rows
  VTM.l.local[is.nan(VTM.l.local)] <- 0L
  VTM.u.local[is.nan(VTM.u.local)] <- 0L

# global
  total.rows <- matrix(rep(colSums(x), K), J)
  muestra_vjk_global <- apply(muestra_vjk_local, c(1L, 2L, 4L), sum)
  VTM.lower <- apply(muestra_vjk_global, c(1L, 2L),
                     stats::quantile, probs = (1L - confidence)/2L)
  VTM.upper <- apply(muestra_vjk_global, c(1L, 2L),
                     stats::quantile, probs = 1L - (1L - confidence)/2L)
  VTM.mean <- apply(muestra_vjk_global, c(1L, 2L), mean)
  inf <- VTM.mean - VTM.lower
  sup <- VTM.upper - VTM.mean
  VTM.votes <- apply(vjk.units, c(1L, 2L), sum)
  VTM.lower <- pmax(pmin(VTM.votes - inf, VTM.votes), 0L)
  VTM.upper <- pmin(pmax(VTM.votes + sup, VTM.votes), total.rows)
  VTM.lower <- VTM.lower/rowSums(VTM.votes)
  VTM.upper <- VTM.upper/rowSums(VTM.votes)

  return(list("VTM.l.local" = VTM.l.local, "VTM.u.local" = VTM.u.local,
              "VTM.lower" = VTM.lower, "VTM.upper" = VTM.upper))
}

