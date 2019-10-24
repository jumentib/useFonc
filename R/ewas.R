#' simulator_ewas : function to simulate DNA methylation data for EWAS
#'
#' @param n	: number of individuals
#' @param p	: number of cpg variables
#' @param K	: number of latent factors
#' @param freq : (vector) mean methylation values (if NULL, set randomly)
#' @param prop.causal : proportion of causal variables (probes/loci)
#' @param prop.variance : proportion of exposure variance explained by latent structure (intensity of confounding)
#' @param sigma :	standard deviation of residual errors
#' @param sd.A :	standard deviation for effect sizes
#' @param mean.A :	(vector) mean of effect sizes
#' @param sd.U : (vector) standard deviations for factors
#' @param sd.V : standard deviations for loadings
#'
#' @return Y : matrix of methylation beta values
#' @return Ynorm : pnorm(Y)
#' @return X : exposure
#' @return A : effect sizes exposure
#' @return causal : set of CpGs associated with the exposure
#' @return U : simulated confounders
#' @return V : loadings of coufounders
#' @return freq : mean methylation values
#'
#' @details
#'
#' This function is used to simulate datasets for EWAS.
#' The simulation model is based on linear relationships.
#' First, it construct a covariance matrix for X and U and prop.variance
#' (intensity of the confounders or correlation between X and U).
#' Then this matrix is used to simulate via normal laws X and U.
#' Thereafter, the effect sizes of X (A) and U (V) are calculated
#' using mean parameters of effect sizes (meanA) and standard deviations (sdA and sdV).
#' Note that the effect sizes of X are calculated only for causal mediators with X.
#' For non-causal mediators, the effect sizes is 0.
#' On the other hand, a residual error matrix is calculated via the sigma (Z) parameter.
#' To finish the methylation matrix is calculated thanks to the formula : Y = V*U + A*X + Z
#' @examples
#' # Simulate data :
#' simu <- simulator_ewas(100, 500, 5)
#' @export
simulator_ewas <- function (n = 100,
                            p = 500,
                            K = 5,
                            freq = NULL,
                            prop.causal = 0.010,
                            prop.variance = 0.4,
                            sigma = 1,
                            sd.A = 0.2,
                            mean.A = 1.0,
                            sd.U = 1.0,
                            sd.V = 1.0)
{
  causal <- sample.int(p, prop.causal * p)

  x.nb = length(causal)

  if (is.null(freq)) freq <- runif(n = p,min =  0.2,max =  0.8) # mean of methylation for each site

  cs <- runif(K, min = -1, max = 1)

  theta <- sqrt( prop.variance /sum((cs/sd.U)^2) )

  # constructing the covariance matrix
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K)

  Sigma <- rbind(Sigma, matrix(cs*theta, nrow = 1))

  Sigma <- cbind(Sigma, matrix(c(cs*theta, 1), ncol = 1))

  UX <- MASS::mvrnorm(n, mu = rep(0, K + 1), Sigma = Sigma)
  U <- UX[, 1:K, drop = FALSE]   # confounders
  X <- UX[, K + 1, drop = FALSE] # outcome

  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))

  A <- matrix(0, p, 1)
  A[causal, 1] <- rnorm(x.nb, mean.A, sd.A)

  Epsilon <- apply(matrix(rep(0,p),nrow = 1), 2, function(x) rnorm(n,x,sigma))

  Z = U %*% t(V) + X %*% t(A) + Epsilon

  M = matrix(rep(qnorm(freq),n), nrow = n, byrow = T) + Z

  N = M
  M = pnorm(M)

  return(list(Y = N,
              X = X,
              A = A,
              U = U,
              V = V,
              Z = Z,
              K = K,
              Ynorm = M,
              causal = causal))
}

#' simulator_ewas_real : function to simulate DNA methylation data for EWAS (looking like real data)
#'
#' @param n	: number of individuals
#' @param p	: number of cpg variables
#' @param K	: number of latent factors
#' @param freq : (vector) mean methylation values (if NULL, set randomly)
#' @param prop.causal : proportion of causal variables (probes/loci)
#' @param prop.variance : proportion of exposure variance explained by latent structure (intensity of confounding)
#' @param sigma :	standard deviation of residual errors
#' @param sd.A :	standard deviation for effect sizes
#' @param mean.A :	(vector) mean of effect sizes
#' @param sd.U : (vector) standard deviations for factors
#' @param sd.V : standard deviations for loadings
#'
#' @return Y : matrix of methylation beta values
#' @return Ynorm : pnorm(Y)
#' @return X : exposure
#' @return A : effect sizes exposure
#' @return causal : set of CpGs associated with the exposure
#' @return U : simulated confounders
#' @return V : loadings of coufounders
#' @return freq : mean methylation values
#'
#' @details
#'
#' This function is used to simulate datasets for EWAS.
#' The simulation model is based on linear relationships.
#' First, it construct a covariance matrix for X and U and prop.variance
#' (intensity of the confounders or correlation between X and U).
#' Then this matrix is used to simulate via normal laws X and U.
#' Thereafter, the effect sizes of X (A) and U (V) are calculated
#' using mean parameters of effect sizes (meanA) and standard deviations (sdA and sdV).
#' Note that the effect sizes of X are calculated only for causal mediators with X.
#' For non-causal mediators, the effect sizes is 0.
#' On the other hand, a residual error matrix is calculated via the sigma (Z) parameter.
#' To finish the methylation matrix is calculated thanks to the formula : Y = V*U + A*X + Z
#' NEED IMPROVEMENT
#'
#' @examples
#' # Simulate data :
#' simu <- simulator_ewas_real(100, 500, 5)
#' @export
simulator_ewas_real <- function (n = 100,
                                 p = 500,
                                 K = 5,
                                 freq = NULL,
                                 sd.freq = 1,
                                 prop.causal = 0.010,
                                 prop.variance = 0.4,
                                 sigma = 1,
                                 sd.A = 0.2,
                                 mean.A = 1.0,
                                 sd.U = 1.0,
                                 sd.V = 1.0,
                                 pour = prop.causal*2)
{
  if (length(sd.freq) == 1) {
    causal <- sample.int(p, prop.causal * p)
  }
  else {
    ord <- order(sd.freq, decreasing = T)[1:(pour*p)] # 10% CpGs most variant

    causal <- sample(ord, prop.causal * p)
  }


  x.nb = length(causal)

  if (is.null(freq)) freq <- runif(n = p,min =  0.2,max =  0.8) # mean of methylation for each site
  #if (is.null(var)) var <- rep(1, p) # mean of methylation for each site

  cs <- runif(K, min = -1, max = 1)

  theta <- sqrt( prop.variance /sum((cs/sd.U)^2) )

  # constructing the covariance matrix
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K)

  Sigma <- rbind(Sigma, matrix(cs*theta, nrow = 1))

  Sigma <- cbind(Sigma, matrix(c(cs*theta, 1), ncol = 1))

  UX <- MASS::mvrnorm(n, mu = rep(0, K + 1), Sigma = Sigma)
  U <- UX[, 1:K, drop = FALSE]   # confounders
  X <- UX[, K + 1, drop = FALSE] # outcome

  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))

  A <- matrix(0, p, 1)
  A[causal, 1] <- rnorm(x.nb, mean.A, sd.A)

  Epsilon <- apply(matrix(rep(0,p),nrow = 1), 2, function(x) rnorm(n,x,sigma))

  Z = U %*% t(V) + X %*% t(A) + Epsilon

  M = matrix(rep(qnorm(freq),n), nrow = n, byrow = T) + Z * (sd.freq)

  #M <- NULL
  #for (i in 1:p) {
  #  m <- rnorm(n, freq[i], sqrt(var[i]))
  #  M <- cbind(M, m)
  #}


  N = M
  M = pnorm(M)

  return(list(Y = N,
              X = X,
              A = A,
              U = U,
              V = V,
              Z = Z,
              K = K,
              Ynorm = M,
              causal = causal))
}
