#' med.test : The Sobel mediation test
#'
#' @description
#'
#' To compute statistics and p-values for the Sobel test. Results for
#' three versions of "Sobel test" are provided: Sobel test, Aroian test and Goodman test.
#' Function adapt from package "bda". (for covariable)
#'
#' @param mv The mediator variable (aka M)
#' @param iv The independant variable (aka X)
#' @param dv The dependent variable (aka Y)
#' @param cof cofacteur
#' @return pValue et score for Sobel test, Aroian test and Goodman test.
#'
#' @details
#'
#' To test whether a mediator carries the influence on an IV to a DV.
#' Missing values are not allowed.
#'
#' @export
med.test <- function(mv, iv, dv, cof = NULL) {
  nm = length(mv)
  ni = length(iv)
  nd = length(dv)

  if (is.null(cof)) {
    tmp = summary(lm(mv ~ iv))
  } else {
    mod1 <- data.frame(mv, iv, cof)
    tmp = summary(lm(mv ~ iv + ., data = mod1))
  }
  a = tmp$coef[2, 1]
  sa = tmp$coef[2, 2]
  if (is.null(cof)) {
    tmp = summary(lm(dv ~ mv + iv))
  } else {
    mod2 <- data.frame(dv, mv, iv, cof)
    tmp = summary(lm(dv ~ mv + iv + ., data = mod2))
  }
  b = tmp$coef[2, 1]
  sb = tmp$coef[2, 2]
  tmp1 = b^2 * sa^2 + a^2 * sb^2
  tmp2 = sa^2 * sb^2
  zsob = a * b/sqrt(tmp1)
  psob = pnorm(-abs(zsob)) * 2
  zaro = a * b/sqrt(tmp1 + tmp2)
  paro = pnorm(-abs(zaro)) * 2
  if (tmp1 > tmp2) {
    zgm = a * b/sqrt(tmp1 - tmp2)
    pgm = pnorm(-abs(zgm)) * 2
  } else {
    zgm = NA
    pgm = NA
  }
  p.value = c(psob, paro, pgm)
  z.value = c(zsob, zaro, zgm)
  out = data.frame(rbind(z.value, p.value))
  names(out) = c("Sobel", "Aroian", "Goodman")
  return(out)

}

#' sob.parallel : The Sobel mediation test for high dimension
#'
#' @description
#'
#' To compute statistics and p-values for the Sobel test. Results for
#' three versions of "Sobel test" are provided: Sobel test, Aroian test and Goodman test.
#' Function adapt from package "bda". (for covariable)
#'
#' @param X exposure
#' @param Y outcome
#' @param M methylation matrix
#' @param conf latent factors
#' @param nb.core number of core used (parallel)
#' @return pValue et score for Sobel test, Aroian test and Goodman test.
#'
#' @details
#'
#' To test whether a mediator carries the influence on an IV to a DV.
#' Missing values are not allowed.
#'
#' @export
sob.parallel <- function(X, Y, M, conf = NULL, nb.core = 2) {
  p <- ncol(M)
  med <- function(i) {
    # argument
    if (is.null(conf)) {
      res <- med.test(mv = M[, i], iv = X, dv = Y, cof = NULL)
    }
    else {
      res <- med.test(mv = M[, i], iv = X, dv = Y, cof = conf)
    }

    return(data.frame(pv.sobel = res[2, 1], pv.aroian = res[2, 3], pv.goodman = res[2, 3],
                      zv.sobel = res[1, 1], zv.aroian = res[1, 3], zv.goodman = res[1, 3]))
  }
  return(do.call("rbind", mclapply(1:p, med, mc.cores = nb.core)))
}



#' r_mediation : function to simulate DNA methylation data for mediation analyzis
#'
#' @param n	: number of individuals
#' @param p	: number of cpg variables
#' @param K	: number of latent factors
#' @param freq : (vector) mean methylation values (if NULL, set randomly)
#' @param prop.causal : proportion of causal variables (probes/loci)
#' @param prop.causal.x : proportion of causal cpg M -> x
#' @param prop.causal.y : proportion of causal cpg M -> y
#' @param prop.causal.ylx : proportion of causal y in causal x
#' @param prop.variance.y : proportion of phenotypic variance explained by latent structure (intensity of confounding)
#' @param prop.variance.x : proportion of exposure variance explained by latent structure (intensity of confounding)
#' @param rho : correlation outcome/exposure (direct effect)
#' @param sigma :	standard deviation of residual errors
#' @param sd.B : standard deviation for effect sizes (B: M->Y)
#' @param mean.B :	(vector) mean of effect sizes
#' @param sd.A :	standard deviation for effect sizes (A: M->X)
#' @param mean.A :	(vector) mean of effect sizes
#' @param sd.U : (vector) standard deviations for factors
#' @param sd.V : standard deviations for loadings
#'
#' @return M : matrix of methylation beta values
#' @return Y : phenotype/health outcome
#' @return B : effect sizes phenotype/health outcome
#' @return X : exposure
#' @return A : effect sizes exposure
#' @return mediators : set of true mediators
#' @return causal.x : set of CpGs associated with the exposure
#' @return causal.y : set of CpGs associated with the outcome
#' @return U : simulated confounders
#' @return V : loadings of coufounders
#' @return freq : mean methylation values
#' @return controls : true control gene (NOT USE for simulation study)
#'
#' @details
#'
#' This function is used to simulate datasets for analysis of mediations.
#' The simulation model is based on linear relationships.
#' First, it construct a covariance matrix for X, Y and U using the parameter rho
#' (direct effect or correlation between X and Y) and propvar
#' (intensity of the confounders or correlation between Y and U).
#' Then this matrix is used to simulate via normal laws X, Y and U.
#' Thereafter, the effect sizes of X (A), Y (B) and U (V) are calculated
#' using mean parameters of effect sizes (meanA and meanB) and standard deviations (sdA, sdB and sdV).
#' Note that the effect sizes of X and Y are calculated only for causal mediators with X and/or Y.
#' For non-causal mediators, the effect sizes is 0.
#' On the other hand, a residual error matrix is calculated via the sigma (Z) parameter.
#' To finish the methylation matrix is calculated thanks to the formula : M = V*U + A*X + B*Y + Z
#' @examples
#' # Simulate data :
#' simu <- r_mediation(100, 500, 5)
#' @export
r_mediation <- function(n,
                        p,
                        K,
                        freq = NULL,
                        prop.causal.x = 0.010,
                        prop.causal.y = 0.010,
                        prop.causal.ylx = 0.5,
                        prop.variance.y = 0.6,
                        prop.variance.x = 0.2,
                        rho = 0.2,
                        sigma = 0.2,
                        sd.A = 1.0,
                        mean.A = 3.0,
                        sd.B = 1.0,
                        mean.B = 5.0,
                        sd.U = 1.0,
                        sd.V = 1.0)
{
  causal.x <- sample.int(p, prop.causal.x * p)
  causal.ylx <- sample(causal.x , prop.causal.ylx*length(causal.x))
  if (prop.causal.y * p < prop.causal.ylx*length(causal.x)) {
    stop("# causal y < # mediators")
  }
  else {
    causal.y <- c(causal.ylx, sample.int(p, prop.causal.y * p - prop.causal.ylx*length(causal.x)) )
  }
  x.nb = length(causal.x)
  y.nb = length(causal.y)

  if (is.null(freq)) freq <- runif(n = p,min =  0.2,max =  0.8) # mean of methylation for each site

  if (prop.variance.y + rho^2 > 1) stop("prop.variance.y + rho^2 > 1")
  if (prop.variance.x + rho^2 > 1) stop("prop.variance.x + rho^2 > 1")

  #cs <- runif(K, min = -1, max = 1)
  #theta.y <- sqrt( prop.variance.y /sum((cs/sd.U)^2) )
  #theta.x <- sqrt( prop.variance.x /sum((cs/sd.U)^2) )

  cs.y <- runif(K, min = -1, max = 1)
  cs.x <- runif(K, min = -1, max = 1)
  theta.y <- sqrt( prop.variance.y /sum((cs.y/sd.U)^2) )
  theta.x <- sqrt( prop.variance.x /sum((cs.x/sd.U)^2) )

  # constructing the covariance matrix
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K)

  Sigma <- rbind(Sigma, matrix(cs.y*theta.y, nrow = 1))
  Sigma <- rbind(Sigma, matrix(cs.x*theta.x, nrow = 1))

  Sigma <- cbind(Sigma, matrix(c(cs.y*theta.y, 1, rho), ncol = 1))
  Sigma <- cbind(Sigma, matrix(c(cs.x*theta.x, rho, 1), ncol = 1))

  UYX <- MASS::mvrnorm(n, mu = rep(0, K + 2), Sigma = Sigma)
  U <- UYX[, 1:K, drop = FALSE]   # confounders
  Y <- UYX[, K + 1, drop = FALSE] # outcome
  X <- UYX[, K + 2, drop = FALSE] # exposure

  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))

  A <- matrix(0, p, 1)
  A[causal.x, 1] <- rnorm(x.nb, mean.A, sd.A)

  B <- matrix(0, p, 1)
  B[causal.y, 1] <- rnorm(y.nb, mean.B, sd.B)

  Epsilon <- apply(matrix(rep(0,p),nrow = 1), 2, function(x) rnorm(n,x,sigma))

  Z = U %*% t(V) + X %*% t(A) + Y %*% t(B) + Epsilon

  M = matrix(rep(qnorm(freq),n), nrow = n, byrow = T) + Z

  M = pnorm(M)

  return(list(M = M,
              Y = Y,
              B = B,
              X = X,
              A = A,
              mediators = sort(causal.ylx),
              causal.x = sort(causal.x),
              causal.y = sort(causal.y),
              U = U,
              V = V,
              freq = freq,
              Sigma = Sigma,
              controls = !(1:p %in% unique(sort(c(sort(causal.x),sort(causal.y)))))))
}



#' r_mediation_real_gene : function to simulate DNA methylation data for mediation analyzis
#'
#' @param meth : methylation matrix
#' @param n	: number of individuals
#' @param p	: number of cpg variables
#' @param K	: number of latent factors
#' @param prop.causal : proportion of causal variables (probes/loci)
#' @param prop.causal.x : proportion of causal cpg M -> x
#' @param prop.causal.y : proportion of causal cpg M -> y
#' @param prop.causal.ylx : proportion of causal y in causal x
#' @param prop.variance.y : proportion of phenotypic variance explained by latent structure (intensity of confounding)
#' @param prop.variance.x : proportion of exposure variance explained by latent structure (intensity of confounding)
#' @param rho : correlation outcome/exposure (direct effect)
#' @param sigma :	standard deviation of residual errors
#' @param sd.B : standard deviation for effect sizes (B: M->Y)
#' @param mean.B :	(vector) mean of effect sizes
#' @param sd.A :	standard deviation for effect sizes (A: M->X)
#' @param mean.A :	(vector) mean of effect sizes
#' @param sd.U : (vector) standard deviations for factors
#' @param sd.V : standard deviations for loadings
#' @param s.real : strengh of real data
#'
#' @return M : matrix of methylation beta values
#' @return Y : phenotype/health outcome
#' @return B : effect sizes phenotype/health outcome
#' @return X : exposure
#' @return A : effect sizes exposure
#' @return mediators : set of true mediators
#' @return causal.x : set of CpGs associated with the exposure
#' @return causal.y : set of CpGs associated with the outcome
#' @return U : simulated confounders
#' @return V : loadings of coufounders
#' @return freq : mean methylation values
#' @return controls : true control gene (NOT USE for simulation study)
#'
#' @details
#'
#' This function is used to simulate datasets for analysis of mediations.
#' The simulation model is based on linear relationships.
#' First, it construct a covariance matrix for X, Y and U using the parameter rho
#' (direct effect or correlation between X and Y) and propvar
#' (intensity of the confounders or correlation between Y and U).
#' Then this matrix is used to simulate via normal laws X, Y and U.
#' Thereafter, the effect sizes of X (A), Y (B) and U (V) are calculated
#' using mean parameters of effect sizes (meanA and meanB) and standard deviations (sdA, sdB and sdV).
#' Note that the effect sizes of X and Y are calculated only for causal mediators with X and/or Y.
#' For non-causal mediators, the effect sizes is 0.
#' On the other hand, a residual error matrix is calculated via the sigma (Z) parameter.
#' To finish the methylation matrix is calculated thanks to the formula : M = meth + V*U + A*X + B*Y + Z
#' @examples
#' # Simulate data :
#'
#' @export
r_mediation_real_gene <- function(meth = NULL, n = nrow(meth), p = ncol(meth), K = 5,
                                  prop.causal.x = 0.005, prop.causal.y = 0.005, prop.causal.ylx = 0.5,
                                  prop.variance.y = 0.1, prop.variance.x = 0.1, rho = 0.2, sigma = 1,
                                  sd.A = 0.1, mean.A = 1, sd.B = 0.1, mean.B = 1,
                                  sd.U = sort(runif(K, 0, 2)), sd.V = 5, s.real = 1) {

  causal.x <- sample.int(p, prop.causal.x * p)
  causal.ylx <- sample(causal.x , prop.causal.ylx*length(causal.x))
  if (prop.causal.y * p < prop.causal.ylx*length(causal.x)) {
    stop("# causal y < # mediators")
  } else {
    causal.y <- c(causal.ylx, sample.int(p, prop.causal.y * p - prop.causal.ylx*length(causal.x)) )
  }
  x.nb = length(causal.x)
  y.nb = length(causal.y)

  # if (is.null(freq)) freq <- runif(n = p,min =  0.2,max =  0.8) # mean of methylation for each site

  if (prop.variance.y + rho^2 > 1) stop("prop.variance.y + rho^2 > 1")
  if (prop.variance.x + rho^2 > 1) stop("prop.variance.x + rho^2 > 1")

  #cs <- runif(K, min = -1, max = 1)
  #theta.y <- sqrt( prop.variance.y /sum((cs/sd.U)^2) )
  #theta.x <- sqrt( prop.variance.x /sum((cs/sd.U)^2) )

  cs.y <- runif(K, min = -1, max = 1)
  cs.x <- runif(K, min = -1, max = 1)
  theta.y <- sqrt(prop.variance.y / sum((cs.y / sd.U)^2) )
  theta.x <- sqrt(prop.variance.x / sum((cs.x / sd.U)^2) )

  # constructing the covariance matrix
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K)

  Sigma <- rbind(Sigma, matrix(cs.y*theta.y, nrow = 1))
  Sigma <- rbind(Sigma, matrix(cs.x*theta.x, nrow = 1))

  Sigma <- cbind(Sigma, matrix(c(cs.y*theta.y, 1, rho), ncol = 1))
  Sigma <- cbind(Sigma, matrix(c(cs.x*theta.x, rho, 1), ncol = 1))

  UYX <- MASS::mvrnorm(n, mu = rep(0, K + 2), Sigma = Sigma)
  U <- UYX[, 1:K, drop = FALSE]   # confounders
  Y <- UYX[, K + 1, drop = FALSE] # outcome
  X <- UYX[, K + 2, drop = FALSE] # exposure

  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))

  A <- matrix(0, p, 1)
  A[causal.x, 1] <- rnorm(x.nb, mean.A, sd.A)

  B <- matrix(0, p, 1)
  B[causal.y, 1] <- rnorm(y.nb, mean.B, sd.B)

  Epsilon <- apply(matrix(rep(0, p), nrow = 1), 2, function(x) rnorm(n, x, sigma))

  Z <- U %*% t(V) + X %*% t(A) + Y %*% t(B) + Epsilon + qnorm(meth) * s.real

  M <- pnorm(Z)

  return(list(M = M,
              X = X,
              Y = Y,
              A = A,
              B = B,
              mediators = sort(causal.ylx),
              causal.x = sort(causal.x),
              causal.y = sort(causal.y),
              U = U,
              V = V,
              Sigma = Sigma))
}


#' r_mediation_natural : function to simulate DNA methylation data for mediation analyzis (based on NaturalGWAS model)
#'
#' @param meth : methylation matrix
#' @param eff.size.x : effect size of X on M
#' @param eff.size.x : effect size of Y on M
#' @param n.x : number of causal CpGs associated with X
#' @param n.y : number of causal CpGs associated with y
#' @param n.xy : number of true mediators
#' @param k : number of latent factors (need to be estimated with a pca)
#' @param chr : chromosome
#' @param pos : chromosome position
#' @param window : distance between causal CpGs
#' @param pca : a former pca with svd() function, speed the simulation (NULL by default)
#'
#' @return causal.x : set of CpGs associated with the exposure
#' @return causal.y : set of CpGs associated with the outcome
#' @return causal.xy : set of true mediators
#' @return X : simulate phenotype
#' @return Y : simulate outcome
#'
#' @details
#'
#' Adaptation of the "simu_pheno" function of the "naturalgwas" package to simulate mediation data in high dimension
#' @examples
#' # Simulate data :
#'
#' @export
r_mediation_natural <- function(meth = NULL, eff.size.x = 100, eff.size.y = 100,
                                n.x = 10, n.y = 10, n.xy = 5, k = NULL,
                                chr = NULL, pos = NULL, window = 101, pca = NULL) {
  require(naturalgwas)

  # pca
  if (is.null(pca)) {
    pc1 <- svd(meth)
  } else {
    pc1 <- pca
  }

  pc1.sdev2 <- pc1$d^2
  sigma <- sqrt(sum(pc1.sdev2) - sum(pc1.sdev2[1:7]))
  base.effect <- sqrt(sum(pc1.sdev2))

  # ref set
  chrpos <- cbind(chr = chr, pos = pos)

  ref.set <- create_refset(chrpos, window = window)

  causal.x <- sort(sample(ref.set, n.x, rep = FALSE))
  #causal.x <- sample.int(p, prop.causal.x * p)
  causal.ylx <- sort(sample(causal.x , n.xy))
  causal.y <- sort(c(causal.ylx, sample(ref.set, n.y - n.xy, rep = FALSE)))

  sigma <- log10(sigma)
  base.effect <- log10(base.effect)

  effect.size.x <- eff.size.x * base.effect
  effect.size.y <- eff.size.y * base.effect

  X <- effect.size.x * rowSums(meth[, causal.x]) +
    rowSums(pc1$u[, 1:k]) + rnorm(nrow(meth), sd = sigma)

  Y <- effect.size.y * rowSums(meth[, causal.y]) +
    rowSums(pc1$u[, 1:k]) + rnorm(nrow(meth), sd = sigma)

  return(list(X = X,
              Y = Y,
              causal.xy = causal.ylx,
              causal.x = causal.x,
              causal.y = causal.y))
}


#' r_mediation_cell_type : function to simulate DNA methylation data for mediation analyzis (and cell type)
#'
#' @param n	: number of individuals
#' @param p	: number of cpg variables
#' @param K	: number of latent factors
#' @param freq : (vector) mean methylation values (if NULL, set randomly)
#' @param prop.causal : proportion of causal variables (probes/loci)
#' @param prop.causal.x : proportion of causal cpg M -> x
#' @param prop.causal.y : proportion of causal cpg M -> y
#' @param prop.causal.ylx : proportion of causal y in causal x
#' @param prop.variance.y : proportion of phenotypic variance explained by latent structure (intensity of confounding)
#' @param prop.variance.x : proportion of exposure variance explained by latent structure (intensity of confounding)
#' @param rho : correlation outcome/exposure (direct effect)
#' @param sigma :	standard deviation of residual errors
#' @param sd.B : standard deviation for effect sizes (B: M->Y)
#' @param mean.B :	(vector) mean of effect sizes
#' @param sd.A :	standard deviation for effect sizes (A: M->X)
#' @param mean.A :	(vector) mean of effect sizes
#' @param sd.U : (vector) standard deviations for factors
#' @param sd.V : standard deviations for loadings
#' @param K.ct : number of cell type
#' @param sd.ct : standard deviations for loadings (cell type)
#' @param alpha : parameter for the dirichlet distribution (for cell type), default : runif(K.ct)
#'
#' @return M : matrix of methylation beta values
#' @return Y : phenotype/health outcome
#' @return B : effect sizes phenotype/health outcome
#' @return X : exposure
#' @return A : effect sizes exposure
#' @return mediators : set of true mediators
#' @return causal.x : set of CpGs associated with the exposure
#' @return causal.y : set of CpGs associated with the outcome
#' @return U : simulated confounders
#' @return V : loadings of coufounders
#' @return freq : mean methylation values
#' @return controls : true control gene (NOT USE for simulation study)
#' @return cell_type : proportion of cell type
#' @return  tcell_type : loadings of cell type
#'
#' @details
#'
#' This function is used to simulate datasets for analysis of mediations.
#' The simulation model is based on linear relationships.
#' First, it construct a covariance matrix for X, Y and U using the parameter rho
#' (direct effect or correlation between X and Y) and propvar
#' (intensity of the confounders or correlation between Y and U).
#' Then this matrix is used to simulate via normal laws X, Y and U.
#' Thereafter, the effect sizes of X (A), Y (B) and U (V) are calculated
#' using mean parameters of effect sizes (meanA and meanB) and standard deviations (sdA, sdB and sdV).
#' Note that the effect sizes of X and Y are calculated only for causal mediators with X and/or Y.
#' For non-causal mediators, the effect sizes is 0.
#' On the other hand, a residual error matrix is calculated via the sigma (Z) parameter.
#' Cell type is simulate with dirichlet distribution
#' To finish the methylation matrix is calculated thanks to the formula : M = V*U + A*X + B*Y + Z
#' @examples
#' # Simulate data :
#' simu <- r_mediation_cell_type(100, 500, 2, 5)
#' @export
r_mediation_cell_type <- function(n,
                                  p,
                                  K,
                                  K.ct,
                                  freq = NULL,
                                  prop.causal.x = 0.010,
                                  prop.causal.y = 0.010,
                                  prop.causal.ylx = 0.5,
                                  prop.variance.y = 0.6,
                                  prop.variance.x = 0.2,
                                  rho = 0.2,
                                  sigma = 0.2,
                                  sd.A = 1.0,
                                  mean.A = 3.0,
                                  sd.B = 1.0,
                                  mean.B = 5.0,
                                  sd.U = 1.0,
                                  sd.V = 1.0,
                                  sd.ct = 1.0,
                                  alpha = NULL)
{
  causal.x <- sample.int(p, prop.causal.x * p)
  causal.ylx <- sample(causal.x , prop.causal.ylx*length(causal.x))
  if (prop.causal.y * p < prop.causal.ylx*length(causal.x)) {
    stop("# causal y < # mediators")
  }
  else {
    causal.y <- c(causal.ylx, sample.int(p, prop.causal.y * p - prop.causal.ylx*length(causal.x)) )
  }
  x.nb = length(causal.x)
  y.nb = length(causal.y)

  if (is.null(freq)) freq <- runif(n = p,min =  0.2,max =  0.8) # mean of methylation for each site

  if (prop.variance.y + rho^2 > 1) stop("prop.variance.y + rho^2 > 1")
  if (prop.variance.x + rho^2 > 1) stop("prop.variance.x + rho^2 > 1")

  #cs <- runif(K, min = -1, max = 1)
  #theta.y <- sqrt( prop.variance.y /sum((cs/sd.U)^2) )
  #theta.x <- sqrt( prop.variance.x /sum((cs/sd.U)^2) )

  cs.y <- runif(K, min = -1, max = 1)
  cs.x <- runif(K, min = -1, max = 1)
  theta.y <- sqrt( prop.variance.y /sum((cs.y/sd.U)^2) )
  theta.x <- sqrt( prop.variance.x /sum((cs.x/sd.U)^2) )

  # constructing the covariance matrix
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K)

  Sigma <- rbind(Sigma, matrix(cs.y*theta.y, nrow = 1))
  Sigma <- rbind(Sigma, matrix(cs.x*theta.x, nrow = 1))

  Sigma <- cbind(Sigma, matrix(c(cs.y*theta.y, 1, rho), ncol = 1))
  Sigma <- cbind(Sigma, matrix(c(cs.x*theta.x, rho, 1), ncol = 1))

  UYX <- MASS::mvrnorm(n, mu = rep(0, K + 2), Sigma = Sigma)
  U <- UYX[, 1:K, drop = FALSE]   # confounders
  Y <- UYX[, K + 1, drop = FALSE] # outcome
  X <- UYX[, K + 2, drop = FALSE] # exposure

  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))

  A <- matrix(0, p, 1)
  A[causal.x, 1] <- rnorm(x.nb, mean.A, sd.A)

  B <- matrix(0, p, 1)
  B[causal.y, 1] <- rnorm(y.nb, mean.B, sd.B)

  Epsilon <- apply(matrix(rep(0,p),nrow = 1), 2, function(x) rnorm(n,x,sigma))

  # rajout cell type
  if (is.null(alpha)) {
    alpha <- runif(K.ct)
  }

  a <- qnorm(gtools::rdirichlet(n = n, alpha = alpha))
  ta <- Rfast::rmvnorm(p, mu = rep(0, K.ct), sigma = sd.ct^2 * diag(K.ct))


  Z = U %*% t(V) + a %*% t(ta) + X %*% t(A) + Y %*% t(B) + Epsilon

  M = matrix(rep(qnorm(freq),n), nrow = n, byrow = T) + Z

  M = pnorm(M)

  return(list(M = M,
              Y = Y,
              B = B,
              X = X,
              A = A,
              cell_type = a,
              tcell_type = ta,
              mediators = sort(causal.ylx),
              causal.x = sort(causal.x),
              causal.y = sort(causal.y),
              U = U,
              V = V,
              freq = freq,
              Sigma = Sigma,
              controls = !(1:p %in% unique(sort(c(sort(causal.x),sort(causal.y)))))))
}

#' bin : Transformation of a continuous data vector into a binary vector
#'
#' @param x	: continuous vector
#' @param prop : probability of success on each trial
#'
#' @return binairy vector
#'
#' @export
bin <- function(x, prop = 0.5) {
  ord <- order(x, decreasing = T)

  s1 <- round(length(x) * prop)

  x[ord[1:s1]] <- 1
  x[ord[(s1 + 1):length(x)]] <- 0

  return(x)
}


#' r_mediation_real_cell_type : function to simulate DNA methylation data for mediation analyzis (and cell type)
#'
#' @param CT : Cell type proportion. null by default. You can use real data (n * K.ct matrix)
#' @param CT.l : Cell type loading. null by default. You can use real data (K.ct * p matrix)
#' @param n	: number of individuals
#' @param p	: number of cpg variables
#' @param K	: number of latent factors
#' @param K.ct : number of cell type
#' @param freq : (vector) mean methylation values (if NULL, set randomly)
#' @param prop.causal : proportion of causal variables (probes/loci)
#' @param prop.causal.x : proportion of causal cpg M -> x
#' @param prop.causal.y : proportion of causal cpg M -> y
#' @param prop.causal.ylx : proportion of causal y in causal x
#' @param prop.variance.y : proportion of phenotypic variance explained by latent structure (intensity of confounding)
#' @param prop.variance.x : proportion of exposure variance explained by latent structure (intensity of confounding)
#' @param rho : correlation outcome/exposure (direct effect)
#' @param sigma :	standard deviation of residual errors
#' @param sd.B : standard deviation for effect sizes (B: M->Y)
#' @param mean.B :	(vector) mean of effect sizes
#' @param sd.A :	standard deviation for effect sizes (A: M->X)
#' @param mean.A :	(vector) mean of effect sizes
#' @param sd.U : (vector) standard deviations for factors
#' @param sd.V : standard deviations for loadings
#' @param sd.ct : standard deviations for loadings (cell type)
#' @param alpha : parameter for the dirichlet distribution (for cell type), default : runif(K.ct)
#' @param prob.bin : if you use binairy exposure (X), probability of success on each trial
#' @param strength : if you use real data for loading of cell type, strength of the real data
#'
#' @return M : matrix of methylation beta values
#' @return X : exposure
#' @return Y : phenotype/health outcome
#' @return A : effect sizes exposure
#' @return B : effect sizes phenotype/health outcome
#' @return M.bin : matrix of methylation beta values, use if you use the binairy exposure
#' @return X.bin : binairy exposure
#' @return CT : proportion of cell type
#' @return CT.l : loading of cell type
#' @return mediators : set of true mediators
#' @return causal.x : set of CpGs associated with the exposure
#' @return causal.y : set of CpGs associated with the outcome
#' @return U : simulated confounders
#' @return V : loadings of coufounders
#' @return freq : mean methylation values
#' @return controls : true control gene (NOT USE for simulation study)
#'
#' @details
#'
#' This function is used to simulate datasets for analysis of mediations.
#' The simulation model is based on linear relationships.
#' First, it construct a covariance matrix for X, Y and U using the parameter rho
#' (direct effect or correlation between X and Y) and propvar
#' (intensity of the confounders or correlation between Y and U).
#' Then this matrix is used to simulate via normal laws X, Y and U.
#' Thereafter, the effect sizes of X (A), Y (B) and U (V) are calculated
#' using mean parameters of effect sizes (meanA and meanB) and standard deviations (sdA, sdB and sdV).
#' Note that the effect sizes of X and Y are calculated only for causal mediators with X and/or Y.
#' For non-causal mediators, the effect sizes is 0.
#' On the other hand, a residual error matrix is calculated via the sigma (Z) parameter.
#' For Cell type (CT) and loading of cell type (CT.l) you can use real data.
#' If Cell type is simulate :
#' Cell type is simulate with dirichlet distribution, loading of cell type simulate via normal laws.
#' To finish the methylation matrix is calculated thanks to the formula :
#' M = VU + CT.l*CT + AX + BY + Z
#' @examples
#' # Simulate data :
#' simu <- r_mediation_cell_type(100, 500, 2, 5)
#' @export
r_mediation_real_cell_type <- function(CT = NULL,
                                       CT.l = NULL,
                                       n = nrow(CT),
                                       p = ncol(CT.l),
                                       K = 2,
                                       K.ct = 5,
                                       freq = NULL,
                                       prop.causal.x = 0.010,
                                       prop.causal.y = 0.010,
                                       prop.causal.ylx = 0.5,
                                       prop.variance.y = 0.1,
                                       prop.variance.x = 0.1,
                                       rho = 0.1,
                                       sigma = 1,
                                       sd.A = 0.1,
                                       mean.A = 1,
                                       sd.B = 0.1,
                                       mean.B = 1,
                                       sd.U = 1.0,
                                       sd.V = 1.0,
                                       strength = 1,
                                       sd.ct = 1,
                                       alpha = NULL,
                                       prob.bin = NULL)
{
  causal.x <- sample.int(p, prop.causal.x * p)
  causal.ylx <- sample(causal.x , prop.causal.ylx*length(causal.x))
  if (prop.causal.y * p < prop.causal.ylx*length(causal.x)) {
    stop("# causal y < # mediators")
  }
  else {
    causal.y <- c(causal.ylx, sample.int(p, prop.causal.y * p - prop.causal.ylx*length(causal.x)) )
  }
  x.nb = length(causal.x)
  y.nb = length(causal.y)

  if (is.null(freq)) freq <- runif(n = p,min =  0.2,max =  0.8) # mean of methylation for each site

  if (prop.variance.y + rho^2 > 1) stop("prop.variance.y + rho^2 > 1")
  if (prop.variance.x + rho^2 > 1) stop("prop.variance.x + rho^2 > 1")

  #cs <- runif(K, min = -1, max = 1)
  #theta.y <- sqrt( prop.variance.y /sum((cs/sd.U)^2) )
  #theta.x <- sqrt( prop.variance.x /sum((cs/sd.U)^2) )

  cs.y <- runif(K, min = -1, max = 1)
  cs.x <- runif(K, min = -1, max = 1)
  theta.y <- sqrt( prop.variance.y /sum((cs.y/sd.U)^2) )
  theta.x <- sqrt( prop.variance.x /sum((cs.x/sd.U)^2) )

  # constructing the covariance matrix
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K)

  Sigma <- rbind(Sigma, matrix(cs.y*theta.y, nrow = 1))
  Sigma <- rbind(Sigma, matrix(cs.x*theta.x, nrow = 1))

  Sigma <- cbind(Sigma, matrix(c(cs.y*theta.y, 1, rho), ncol = 1))
  Sigma <- cbind(Sigma, matrix(c(cs.x*theta.x, rho, 1), ncol = 1))

  UYX <- MASS::mvrnorm(n, mu = rep(0, K + 2), Sigma = Sigma)
  U <- UYX[, 1:K, drop = FALSE]   # confounders
  Y <- UYX[, K + 1, drop = FALSE] # outcome
  X <- UYX[, K + 2, drop = FALSE] # exposure

  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))

  A <- matrix(0, p, 1)
  A[causal.x, 1] <- rnorm(x.nb, mean.A, sd.A)

  B <- matrix(0, p, 1)
  B[causal.y, 1] <- rnorm(y.nb, mean.B, sd.B)

  Epsilon <- apply(matrix(rep(0,p),nrow = 1), 2, function(x) rnorm(n,x,sigma))


  # rajout cell type

  if (is.null(CT)) {
    if (is.null(alpha)) {
      alpha <- runif(K.ct)
    }
    CT <- qnorm(gtools::rdirichlet(n = n, alpha = alpha))
  } else {
    # No Zero and no 1
    CT[which(CT <= 0, arr.ind = T)] <- 0.00000001
    CT[which(CT >= 1, arr.ind = T)] <- 0.99999999
    CT <- qnorm(CT)
  }

  if (is.null(CT.l)) {
    CT.l <- t(Rfast::rmvnorm(p, mu = rep(0, ncol(CT)), sigma = sd.ct^2 * diag(ncol(CT))))
  } else {
    CT.l[which(CT.l <= 0, arr.ind = T)] <- 0.00000001
    CT.l[which(CT.l >= 1, arr.ind = T)] <- 0.99999999
    CT.l <- qnorm(CT.l) * strength
  }

  # continue
  Z <- U %*% t(V) + CT %*% CT.l + X %*% t(A) + Y %*% t(B) + Epsilon
  M <- matrix(rep(qnorm(freq),n), nrow = n, byrow = T) + Z
  M <- pnorm(M)

  # si binaire
  if (is.null(prob.bin)) {
    X1 <- pbinom(X, 0, 0.5)
  } else {
    X1 <- bin(X, prop = prob.bin)
  }

  X1 <- 2 * X1 - 1

  Z1 <- U %*% t(V) + CT %*% CT.l + X1 %*% t(A) + Y %*% t(B) + Epsilon
  M1 <- matrix(rep(qnorm(freq),n), nrow = n, byrow = T) + Z1
  M1 <- pnorm(M1)


  return(list(M = M,
              X = X,
              Y = Y,
              A = A,
              B = B,
              M.bin = M1,
              X.bin = X1,
              CT = CT,
              CT.l = CT.l,
              mediators = sort(causal.ylx),
              causal.x = sort(causal.x),
              causal.y = sort(causal.y),
              U = U,
              V = V,
              freq = freq,
              Sigma = Sigma,
              controls = !(1:p %in% unique(sort(c(sort(causal.x),sort(causal.y)))))))
}
