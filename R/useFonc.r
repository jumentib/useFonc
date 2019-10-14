#' cor.paired : Correlation calculation with nearest neighbor
#'
#' @param meth matrice of mDNA (the matrix needs to be ordered according to the chromosomes and the chromosomal position)
#' @param pas 1 by default, corresponds to the nearest neighbor. 2 for the second closest neighbor and etc ...
#' @return a vector of correlation
#'
#' @details
#'
#' This function makes it possible to calculate the correlation with the nearest neighbor (or the second or third, etc.). Correlation of pearson by default.
#'
#' @export
cor.paired <- function(meth, pas = 1) {
  c <- NULL
  for (i in 1 : (ncol(meth) - pas)) {
    c1 <- cor(meth[, i], meth[, i + pas])
    c <- c(c, c1)
  }
  return(c)
}


#' cor.pas : Correlation calculation with nearest neighbor (with genetic distance)
#'
#' @param meth matrice of mDNA (the matrix needs to be ordered according to the chromosomes and the chromosomal position)
#' @param pas 1 by default, corresponds to the nearest neighbor. 2 for the second closest neighbor and etc ...
#' @return a vector of correlation
#'
#' @details
#'
#' This function makes it possible to calculate the correlation with the nearest neighbor (or the second or third, etc.). Correlation of pearson by default.
#' The function cor.pas and different from the functions cor.paired, for not > 1 it does not calculate the correlation with all the neighbors (see the code of the function to understand).
#' @export
cor.pas <- function(meth, pas = 1) {
  meth <- meth[, seq(1, ncol(meth), pas)]
  c <- NULL
  for (i in 1 : (ncol(meth) - 1)) {
    c1 <- cor(meth[, i], meth[, i + 1])
    c <- c(c, c1)
  }
  return(c)
}


#' cor.ld : Correlation calculation with nearest neighbor in mDNA matrix
#'
#' @param meth matrice of mDNA
#' @param annot annotation file of the mDNA matrix (Illuminina type)
#' @param pas 1 by default, corresponds to the nearest neighbor. 2 for the second closest neighbor and etc ...
#' @return a matrix, with a correlation vector, a genetic distance vector.
#'
#' @details
#'
#' This function makes it possible to calculate the correlation with the nearest neighbor (or the second or third, etc.). Correlation of pearson by default.
#' This function orders the mDNA matrix following chromosomes and chromosome
#' positions and uses the cor.pas() function to calculate the correlation for each chromosome.
#' Moreover the distance between each neighbor is calculated.
#'
#' @export
cor.ld <- function(meth, annot, pas = 1) {

  sc <- data.frame(colnames(meth))

  annot <- merge.data.frame(sc, annot, by = 1)
  # connversion chr en numeric
  annot$CHR[annot$CHR == "X"] <- 23
  annot$CHR[annot$CHR == "Y"] <- 24


  annot$CHR <- as.numeric(annot$CHR)

  # annotation ranger dans l'ordre
  a <- doBy::orderBy(~ CHR + MAPINFO, annot)
  cpg <- data.frame(id = a$Name, chr = a$CHR, pos = a$MAPINFO,
                    gene = a$UCSC_RefGene_Name, index = 1:nrow(a))
  cpg <- na.omit(cpg)

  # mDNA dans l'ordre
  meth <- meth[,cpg$id]
  #print(dim(meth))

  # liste des chr
  u <- unique(cpg$chr)

  cdc <- NULL
  for (j in u) {

    print(paste0("CHR : ", j))
    # juste le premier chr
    mi <- meth[, cpg$chr == j]

    # correlation dans le premier chr
    ci <- cor.pas(mi, pas)

    # position sur le chr
    posi <- cpg$pos[cpg$chr == j]
    posi <- posi[seq(1, ncol(mi), pas)]

    # calcul distance
    di <- NULL
    for (i in 2:length(posi)) {
      di <- c(di, posi[i] - posi[i-1])
    }

    cdci <- cbind(cor = ci, dist = di, chr = rep(j, length(ci)),
                  pas = rep(pas, length(ci)))
    #print(dim(cdci))
    cdc <- rbind(cdc, cdci)
  }
  return(cdc)
}

#' cor.ld.pas : Correlation calculation with nearest neighbor in mDNA matrix (for multiple pas)
#'
#' @param meth matrice of mDNA
#' @param annot annotation file of the mDNA matrix (Illuminina type)
#' @param pas 1 by default, corresponds to the nearest neighbor. 2 for the second closest neighbor and etc ...
#' @return a matrix, with a correlation vector, a genetic distance vector.
#'
#' @details
#'
#' This function makes it possible to calculate the correlation with the nearest neighbor (or the second or third, etc.). Correlation of pearson by default.
#' This function orders the mDNA matrix following chromosomes and chromosome
#' positions and uses the cor.pas() function to calculate the correlation for each chromosome.
#' Moreover the distance between each neighbor is calculated.
#' The difference with cor.ld() is that cor.ld.pas() can be used by varying the "pas" (between nearest neighbors).
#'
#' @export
cor.ld.pas <- function(meth, annot, mulpas = 1:5) {
  cdc <- NULL
  for (i in mulpas) {
    print(paste0("-------pas = ", i))
    cdc <- rbind(cdc, cor.ld(meth, annot, pas = i))
  }
  return(cdc)
}


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



#' rank.pwer : Calculate the power and the calibration of statiscal test (for simulation)
#'
#' @param pval pValues from a test
#' @param known.mediator true result of the simulation
#' @param toplist number of high pValues to use to calculate rank power
#' @param decreasing default F. This function can be use with effects size and with decreasing = T
#' @param ral default 0.05. Threshold for the calculation of the power and the recall
#' @return cur : the rank power curve
#' @return auc : the AUC of the rank power curve
#' @return auc.max : the maximum AUC for a perfect model
#' @return auc.norm : auc/auc.max. Normalization of the AUC of the rank power curve
#' @return precision : statistical precision of the test : TP / (TP + FP)
#' @return recall : statistical recall of the test : TP / (TP + FN)
#' @return f1_score : F1 Score of the test : 2(Pre * Rec) / (Pre + Rec)
#' @return length.list : number of pValues lower than "ral" (Bonferroni correction).
#' @return TP : True positive
#' @return FP : False positive
#' @return FN : False negative
#'
#' @details
#'
#' AUC corresponds to the air under the curve (AUC) of a hits rank curve.
#' This curve is obtained by looking for hits in a toplist of the 20 CpGs with the lowest pValues.
#' F1 = 2*recall*power/(recall + power)
#' @export
rank.pwer <- function(pval,known.mediator=NULL,toplist=20, decreasing = F, ral = 0.05) {
  # rank power calculation
  if (is.null(known.mediator) == FALSE) {
    top <- order(pval, decreasing = decreasing)[1:toplist]
    cur <- 0 # rank power curve

    for (i in 1:toplist) {
      cur <- c(cur,sum(top[1:i] %in% known.mediator))
    }
    auc <- sum(cur) # AUC of rank power curve
  }
  else {
    cur <- NA
    auc <- NA
  }

  # AUC standardization
  NH <- length(known.mediator) # number of true hits

  auc.max <- sum(c(1:NH,rep(NH,toplist-NH)))

  al <- ral/length(pval)
  pos <- which(pval <= al)
  neg <- which(pval > al)

  TP <- sum(pos %in% causal)
  FP <- sum(!(pos %in% causal))
  FN <- sum(neg %in% causal)

  pre <- TP / (TP + FP)
  rec <- TP / (TP + FN)
  # taille de list
  l.pv <- length(pos)

  return(list(auc  = auc,
              cur  = cur,
              auc.norm = auc/auc.max,
              auc.max  = auc.max,
              precision = pre,
              recall = rec,
              f1_score = 2*rec*pre/(rec + pre),
              length.list = l.pv,
              TP = TP,
              FP = FP,
              FN = FN))
}


# Write the phenotype data to a file in the format used by GEMMA. Each
# line of the file contains one phenotype observation.

#' write.gemma.pheno
#'
#' @param pheno phenotype
#' @return nothing
#'
#' @export
write.gemma.pheno <- function (file, pheno) {
  phenotype <- 1:nrow(pheno)
  y <- pheno[phenotype] # j'ai retiré une "[]"
  if (is.numeric(y))
    y <- round(y,digits = 6)
  write.table(y,file,quote = FALSE,row.names = FALSE,col.names = FALSE)
}


# Write the covariate data to a file in the format used by GEMMA. Each
# line corresponds to a sample. We must include an additional
# covariate for the intercept.

#' write.gemma.covariates
#'
#' @param covariates covariable
#' @return nothing
#'
#' @export
write.gemma.covariates <- function (file, covariates, pheno) {
  if (length(covariates) == 0) {
    write.table(data.frame(rep(1,nrow(pheno))),file,sep = " ",
                quote = FALSE,row.names = FALSE,col.names = FALSE)
  } else {
    round.col <- function (x) {
      if (is.numeric(x))
        round(x,digits = 6)
      else
        x
    }
    write.table(cbind(1,data.frame(lapply(pheno[covariates],round.col))),
                file,sep = " ",quote = FALSE,row.names = FALSE,
                col.names = FALSE)
  }
}

# Write the SNP information to a space-delimited text file in the
# format used by GEMMA. This file contains one line per SNP, with
# three columns: (1) SNP label, (2) base-pair position, (3)
# chromosome.

#' write.gemma.map
#'
#' @param map map
#' @return nothing
#'
#' @export
write.gemma.map <- function (file, map)
  write.table(map,file,sep = " ",quote = FALSE,
              row.names = FALSE,col.names = FALSE)

# Store the mean genotypes as a space-delimited text file in the
# format used by GEMMA, in which we have one row per SNP, and one
# column per sample. The first three column give the SNP label, and
# the two alleles.

#' write.gemma.geno
#'
#' @param geno genotype
#' @return nothing
#'
#' @export
write.gemma.geno <- function (file, geno, map) {
  geno <- t(geno)
  geno <- as.data.frame(geno,check.names = FALSE)
  geno <- round(geno,digits = 3)
  geno <- cbind(map,geno)
  write.table(geno,file,sep = " ",quote = FALSE,row.names = FALSE,
              col.names = FALSE)
}

# covar : covariable et NULL marche
# pheno : phenotype d'interet (doit etre une matrice)
# geno : genotype --> one row per SNP, and one column per sample.
# map : This file contains one line per SNP, with
# three columns: (1) SNP label, (2) base-pair position, (3)
# chromosome.
# gemma.exe : the executable gemma
# opt.mat : option for the calcul of the related.matrix,
# "1" calculates the centered relatedness # matrix
# while "2" calculates the standardized relatedness matrix
# met.reg : linear mixed model of Bayesian sparse linear mixed model
# lmm : if lmm choose option : where the  "1" performs Wald test
# "2" performs likelihood ratio test, "3" performs score test, "4" performs all the three tests.
# bslmm : if bslmm choose option : "1" fits a standard linear BSLMM
# "2" fits a ridge regression/GBLUP, "3" fits a probit BSLMM

#' gemma : GEMMA programme in R
#'
#' @param covar : covariable et NULL marche
#' @param pheno : phenotype d'interet (doit etre une matrice)
#' @param geno : genotype --> one row per SNP, and one column per sample.
#' @param map : This file contains one line per SNP, with
#' @param three columns: (1) SNP label, (2) base-pair position, (3) chromosome.
#' @param gemma.exe : the executable gemma
#' @param opt.mat : option for the calcul of the related.matrix, "1" calculates the centered relatedness # matrix while "2" calculates the standardized relatedness matrix
#' @param met.reg : linear mixed model of Bayesian sparse linear mixed model.
#' @param lmm : if lmm choose option : where the  "1" performs Wald test "2" performs likelihood ratio test, "3" performs score test, "4" performs all the three tests.
#' @param bslmm : if bslmm choose option : "1" fits a standard linear BSLMM "2" fits a ridge regression/GBLUP, "3" fits a probit BSLMM.
#' @return : a gemma result
#'
#' @details
#'
#' See paper of GEMMA
#'
#' @export
gemma <- function(geno, pheno, map, covar = NULL, gemma.exe = "./gemma.macosx",
                  opt.mat = c(1,2), met.reg = c("lmm", "bslmm"),
                  lmm = c(1,2,3,4), bslmm = c(1,2,3)) {

  system("mkdir gemma")

  write.gemma.pheno("gemma/pheno.txt",pheno)
  write.gemma.covariates("gemma/covariates.txt",covar,pheno)
  write.gemma.geno("gemma/geno.txt",geno,map)
  write.gemma.map("gemma/map.txt",map)

  # related.matrix
  system(paste(gemma.exe,"-g gemma/geno.txt -p gemma/pheno.txt -c gemma/covariates.txt -gk",
               opt.mat,"-o kin"),
         ignore.stdout = F)

  # lmm
  if (met.reg == "lmm") {
    print(met.reg)
    system(paste(gemma.exe,"-g gemma/geno.txt -p gemma/pheno.txt -c gemma/covariates.txt",
                 " -a gemma/map.txt -k output/kin.cXX.txt -lmm", lmm,"-o lmm"),
           ignore.stdout = F)
    res <- read.table("output/lmm.assoc.txt", header = T)
  }
  # bslmm
  else {
    print(met.reg)
    system(paste(gemma.exe,"-g gemma/geno.txt -p gemma/pheno.txt",
                 " -a gemma/map.txt -k output/kin.cXX.txt -bslmm", bslmm," -o bslmm"),
           ignore.stdout = F)
    res <- read.table("output/bslmm.param.txt", header = T)
  }
  return(res)
}

#' gemma2 : GEMMA programme in R (better version of gemma)
#'
#' @param covar : covariable et NULL marche
#' @param pheno : phenotype d'interet (doit etre une matrice)
#' @param geno : genotype --> one row per SNP, and one column per sample.
#' @param map : This file contains one line per SNP, with
#' @param three columns: (1) SNP label, (2) base-pair position, (3) chromosome.
#' @param gemma.exe : the executable gemma
#' @param opt.mat : option for the calcul of the related.matrix, "1" calculates the centered relatedness # matrix while "2" calculates the standardized relatedness matrix
#' @param met.reg : linear mixed model of Bayesian sparse linear mixed model.
#' @param lmm : if lmm choose option : where the  "1" performs Wald test "2" performs likelihood ratio test, "3" performs score test, "4" performs all the three tests.
#' @param bslmm : if bslmm choose option : "1" fits a standard linear BSLMM "2" fits a ridge regression/GBLUP, "3" fits a probit BSLMM.
#' @param pmin : for bslmm : initial no zero proportion of beta (pi)
#' @param pmax : for bslmm : max no zero proportion of beta (pi)
#' @param opt.bs : option for bslmm : w and s default as in Zeng 2017 (10000 for both)
#' @param smax : Number of variants with sparse effects (~ number of major effect loci). Default to 300 like in gemma
#' @return : a gemma result
#'
#' @details
#'
#' See paper of GEMMA
#'
#' @export
gemma2 <- function(geno, pheno, map, covar = NULL, gemma.exe = "./gemma.macosx",
                   opt.mat = c(1,2), met.reg = c("lmm", "bslmm"),
                   lmm = c(1,2,3,4), bslmm = c(1,2,3), pmin = 0.01, pmax = 0.05,
                   opt.bs = "-w  10000 -s  10000", smax = 300) {

  system("mkdir gemma")

  write.gemma.pheno("gemma/pheno.txt",pheno)
  write.gemma.covariates("gemma/covariates.txt",covar,pheno)
  write.gemma.geno("gemma/geno.txt",geno,map)
  write.gemma.map("gemma/map.txt",map)

  # related.matrix
  system(paste(gemma.exe,"-g gemma/geno.txt -p gemma/pheno.txt -c gemma/covariates.txt -gk",
               opt.mat,"-o kin"),
         ignore.stdout = F)

  # lmm
  if (met.reg == "lmm") {
    print(met.reg)
    system(paste(gemma.exe,"-g gemma/geno.txt -p gemma/pheno.txt -c gemma/covariates.txt",
                 " -a gemma/map.txt -k output/kin.cXX.txt -lmm", lmm,"-o lmm"),
           ignore.stdout = F)
    res <- read.table("output/lmm.assoc.txt", header = T)
  }
  # bslmm
  else {
    print(met.reg)
    # opt.bs <- "-w  10000 -s  10000" # Zeng 2017
    pmin <- log10(pmin)
    pmax <- log10(pmax)
    system(paste(gemma.exe,"-g gemma/geno.txt -p gemma/pheno.txt",
                 " -a gemma/map.txt -k output/kin.cXX.txt -bslmm", bslmm," -o bslmm",
                 opt.bs, "-pmin", pmin, "-pmax", pmax, "-smax", smax),
           ignore.stdout = F)
    res <- read.table("output/bslmm.param.txt", header = T)
  }
  return(res)
}

#' pwer.effect : Calculate the power and the calibration of statiscal test (on the effects sizes)
#'
#' @param effect effects sizes from a model
#' @param causal true result of the simulation
#' @param threshold the different thresholds where the power and fdr will be calculated
#' @return the power and the fdr for each threshold
#'
#' @export
pwer.effect <- function(effect, causal, threshold = seq(0, 0.9, 0.05)){
  power <- NULL
  fdr <- NULL

  for (i in threshold) {
    ls <- which(effect > i)
    # power
    pow <- sum(ls %in% causal) / length(causal)
    # recall
    rec <- sum(causal %in% ls) / length(ls)

    # data
    power <- c(power, pow)
    fdr <- c(fdr, rec)

  }
  return(data.frame(threshold, power, fdr))
}

#' lfmm_lasso2 : lfmm lasso with no scale on data
#'
#' @export
lfmm_lasso2 <- function(Y, X, K,
                        nozero.prop = 0.01,
                        lambda.num = 100,
                        lambda.min.ratio = 0.01,
                        lambda = NULL,
                        it.max = 100,
                        relative.err.epsilon = 1e-6) {
  ## init
  m <- lfmm::lassoLFMM(K = K,
                       nozero.prop = nozero.prop,
                       lambda.num = lambda.num,
                       lambda.min.ratio = lambda.min.ratio,
                       lambda = lambda)
  # dat <- LfmmDat(Y = scale(Y, scale = FALSE), X = scale(X, scale = FALSE))
  dat <- lfmm::LfmmDat(Y = Y, X = X) # sans scale
  ## run
  m <- lfmm::lfmm_fit(m, dat, it.max = it.max, relative.err.epsilon = relative.err.epsilon)

  ## return
  m
}

#' lfmm_ridge2 : lfmm ridge with no scale on data
#'
#' @export
lfmm_ridge2 <- function(Y, X, K, lambda = 1e-5, algorithm = "analytical",
                        it.max = 100, relative.err.min = 1e-6) {

  ## init
  m <- lfmm::ridgeLFMM(K = K, lambda = lambda)
  m$algorithm <- algorithm[1]
  # dat <- LfmmDat(Y = scale(Y, scale = FALSE), X = scale(X, scale = FALSE))
  dat <- lfmm::LfmmDat(Y = Y, X = X) # sans scale

  ## run
  m <- lfmm::lfmm_fit(m, dat, it.max = it.max, relative.err.min = relative.err.min)

  ## return
  m
}

#' fn : imput mean for NA data
#'
#' @param x a matrix
#'
#' @export
fn <- function(x){

  m. <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- m.
  return(x)
}

#' panel.hist : pairs function
#'
#' @param x a data.frame
#'
#' @export
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

#' panel.cor : pairs function
#'
#' @param x a data.frame
#'
#' @export
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


#' f1.score.prs : Calcul of the F1 score and the Polygenic Risk Score (PRS) for sparse regression
#'
#' @param sparse.effect : effect from a sparse regression
#' @param causal : true causal SNP
#' @param genotype : genotype --> one row per SNP, and one column per sample.
#' @param phenotype : phenotype of interest
#'
#' @return Return the F1 Score, the precision, the recall, the prediction of the phenotype and the r2 between the phenotype and x.pred.
#'
#' @export
f1.score.prs <- function(sparse.effect, causal, genotype, phenotype) {

  nozero.eff <- which(sparse.effect != 0)
  zero.eff   <- which(sparse.effect == 0)
  # F1 score
  TP <- sum(nozero.eff %in% causal)
  FP <- sum(!(nozero.eff %in% causal))
  FN <- sum(causal %in% zero.eff)

  preci <- TP / (TP + FP)
  recal <- TP / (TP + FN)

  F1 <- 2 * (preci * recal) / (preci + recal)

  # PRS
  x.pred <- genotype[, nozero.eff] %*% matrix(sparse.effect[nozero.eff])
  r2 <- summary(lm(x.pred ~ phenotype))$r.squared
  return(list(f1.score = F1, r2 = r2, precision = preci, recall = recal, x.pred = x.pred))
}



#' rank.pwer.sparse : Calcul the AUC of the Hits Rank Curve for sparse regresion
#'
#' @param sparse.effect : effect from a sparse regression
#' @param causal : true causal SNP
#'
#' @return Return the AUC of the Hits Rank Curve
#'
#' @export
rank.pwer.sparse <- function(sparse.effect, causal) {
  # function
  n.nozero <- sum(sparse.effect != 0)
  selec <- order(abs(sparse.effect), decreasing = T)[1:n.nozero]
  cur <- 0 # rank power curve

  for (i in 1:n.nozero) {
    cur <- c(cur,sum(selec[1:i] %in% causal))
  }
  auc <- sum(cur) # AUC of rank power curve

  # calcul AUC max
  nozero.eff <- which(sparse.effect != 0)
  TP <- sum(nozero.eff %in% causal)
  FP <- sum(!(nozero.eff %in% causal))

  cur.max <- c(1:TP, rep(TP, FP))
  auc.max <- sum(cur.max)

  return(list(auc.norm = auc/auc.max, auc = auc, cur = cur))
}


#' cbind.na : cbind.na
#'
#' @export
cbind.na <- function (..., deparse.level = 1)
{
  na <- nargs() - (!missing(deparse.level))
  deparse.level <- as.integer(deparse.level)
  stopifnot(0 <= deparse.level, deparse.level <= 2)
  argl <- list(...)
  while (na > 0 && is.null(argl[[na]])) {
    argl <- argl[-na]
    na <- na - 1
  }
  if (na == 0)
    return(NULL)
  if (na == 1) {
    if (isS4(..1))
      return(cbind2(..1))
    else return(matrix(...))  ##.Internal(cbind(deparse.level, ...)))
  }
  if (deparse.level) {
    symarg <- as.list(sys.call()[-1L])[1L:na]
    Nms <- function(i) {
      if (is.null(r <- names(symarg[i])) || r == "") {
        if (is.symbol(r <- symarg[[i]]) || deparse.level ==
            2)
          deparse(r)
      }
      else r
    }
  }
  ## deactivated, otherwise no fill in with two arguments
  if (na == 0) {
    r <- argl[[2]]
    fix.na <- FALSE
  }
  else {
    nrs <- unname(lapply(argl, nrow))
    iV <- sapply(nrs, is.null)
    fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
    ## deactivated, otherwise data will be recycled
    #if (fix.na) {
    #    nr <- max(if (all(iV)) sapply(argl, length) else unlist(nrs[!iV]))
    #    argl[[na]] <- cbind(rep(argl[[na]], length.out = nr),
    #        deparse.level = 0)
    #}
    if (deparse.level) {
      if (fix.na)
        fix.na <- !is.null(Nna <- Nms(na))
      if (!is.null(nmi <- names(argl)))
        iV <- iV & (nmi == "")
      ii <- if (fix.na)
        2:(na - 1)
      else 2:na
      if (any(iV[ii])) {
        for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i)))
          names(argl)[i] <- nmi
      }
    }

    ## filling with NA's to maximum occuring nrows
    nRow <- as.numeric(sapply(argl, function(x) NROW(x)))
    maxRow <- max(nRow, na.rm = TRUE)
    argl <- lapply(argl, function(x)  if (is.null(nrow(x))) c(x, rep(NA, maxRow - length(x)))
                   else rbind.na(x, matrix(, maxRow - nrow(x), ncol(x))))
    r <- do.call(cbind, c(argl[-1L], list(deparse.level = deparse.level)))
  }
  d2 <- dim(r)
  r <- cbind2(argl[[1]], r)
  if (deparse.level == 0)
    return(r)
  ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
  ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
  if (ism1 && ism2)
    return(r)
  Ncol <- function(x) {
    d <- dim(x)
    if (length(d) == 2L)
      d[2L]
    else as.integer(length(x) > 0L)
  }
  nn1 <- !is.null(N1 <- if ((l1 <- Ncol(..1)) && !ism1) Nms(1))
  nn2 <- !is.null(N2 <- if (na == 2 && Ncol(..2) && !ism2) Nms(2))
  if (nn1 || nn2 || fix.na) {
    if (is.null(colnames(r)))
      colnames(r) <- rep.int("", ncol(r))
    setN <- function(i, nams) colnames(r)[i] <<- if (is.null(nams))
      ""
    else nams
    if (nn1)
      setN(1, N1)
    if (nn2)
      setN(1 + l1, N2)
    if (fix.na)
      setN(ncol(r), Nna)
  }
  r
}


#' rbind.na : rbind.na
#'
#' @export
rbind.na <- function (..., deparse.level = 1)
{
  na <- nargs() - (!missing(deparse.level))
  deparse.level <- as.integer(deparse.level)
  stopifnot(0 <= deparse.level, deparse.level <= 2)
  argl <- list(...)
  while (na > 0 && is.null(argl[[na]])) {
    argl <- argl[-na]
    na <- na - 1
  }
  if (na == 0)
    return(NULL)
  if (na == 1) {
    if (isS4(..1))
      return(rbind2(..1))
    else return(matrix(..., nrow = 1)) ##.Internal(rbind(deparse.level, ...)))
  }
  if (deparse.level) {
    symarg <- as.list(sys.call()[-1L])[1L:na]
    Nms <- function(i) {
      if (is.null(r <- names(symarg[i])) || r == "") {
        if (is.symbol(r <- symarg[[i]]) || deparse.level ==
            2)
          deparse(r)
      }
      else r
    }
  }

  ## deactivated, otherwise no fill in with two arguments
  if (na == 0) {
    r <- argl[[2]]
    fix.na <- FALSE
  }
  else {
    nrs <- unname(lapply(argl, ncol))
    iV <- sapply(nrs, is.null)
    fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
    ## deactivated, otherwise data will be recycled
    #if (fix.na) {
    #    nr <- max(if (all(iV)) sapply(argl, length) else unlist(nrs[!iV]))
    #    argl[[na]] <- rbind(rep(argl[[na]], length.out = nr),
    #        deparse.level = 0)
    #}
    if (deparse.level) {
      if (fix.na)
        fix.na <- !is.null(Nna <- Nms(na))
      if (!is.null(nmi <- names(argl)))
        iV <- iV & (nmi == "")
      ii <- if (fix.na)
        2:(na - 1)
      else 2:na
      if (any(iV[ii])) {
        for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i)))
          names(argl)[i] <- nmi
      }
    }

    ## filling with NA's to maximum occuring ncols
    nCol <- as.numeric(sapply(argl, function(x) if (is.null(ncol(x))) length(x)
                              else ncol(x)))
    maxCol <- max(nCol, na.rm = TRUE)
    argl <- lapply(argl, function(x)  if (is.null(ncol(x))) c(x, rep(NA, maxCol - length(x)))
                   else cbind(x, matrix(, nrow(x), maxCol - ncol(x))))

    ## create a common name vector from the
    ## column names of all 'argl' items
    namesVEC <- rep(NA, maxCol)
    for (i in 1:length(argl)) {
      CN <- colnames(argl[[i]])
      m <- !(CN %in% namesVEC)
      namesVEC[m] <- CN[m]
    }

    ## make all column names from common 'namesVEC'
    for (j in 1:length(argl)) {
      if (!is.null(ncol(argl[[j]]))) colnames(argl[[j]]) <- namesVEC
    }

    r <- do.call(rbind, c(argl[-1L], list(deparse.level = deparse.level)))
  }

  d2 <- dim(r)

  ## make all column names from common 'namesVEC'
  colnames(r) <- colnames(argl[[1]])

  r <- rbind2(argl[[1]], r)

  if (deparse.level == 0)
    return(r)
  ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
  ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
  if (ism1 && ism2)
    return(r)
  Nrow <- function(x) {
    d <- dim(x)
    if (length(d) == 2L)
      d[1L]
    else as.integer(length(x) > 0L)
  }
  nn1 <- !is.null(N1 <- if ((l1 <- Nrow(..1)) && !ism1) Nms(1))
  nn2 <- !is.null(N2 <- if (na == 2 && Nrow(..2) && !ism2) Nms(2))
  if (nn1 || nn2 || fix.na) {
    if (is.null(rownames(r)))
      rownames(r) <- rep.int("", nrow(r))
    setN <- function(i, nams) rownames(r)[i] <<- if (is.null(nams))
      ""
    else nams
    if (nn1)
      setN(1, N1)
    if (nn2)
      setN(1 + l1, N2)
    if (fix.na)
      setN(nrow(r), Nna)
  }
  r
}

#' RMSE : Calcul Root Square Mean Error
#'
#' @param beta : effect from a linear regression
#' @param true : true effect from a simulation
#'
#' @return the RMSE
#'
#' @export
RMSE <- function(beta, true.beta) {
  return(sqrt(mean((true.beta - beta)^2)))
}

#' F1 : Calcul F1 score for a simulation (on effect sizes)
#'
#' @param beta : effect from a regression
#' @param causal : true causal variable (SNP or CpG)
#' @param nb.hit : number of first hit for calculate the F1 score
#'
#' @return the F1 score, the precision and the recall
#'
#' @export
F1 <- function(beta, causal, nb.hit = 10) {
  abs.beta <- abs(beta)
  p <- length(beta)
  # create list for F1
  first.hit <- order(abs.beta, decreasing = T)[1:nb.hit]
  last.hit <- order(abs.beta, decreasing = T)[(nb.hit + 1):p]
  # F1 score
  tp <- sum(first.hit %in% causal)
  fp <- sum(!(first.hit %in% causal))
  fn <- sum(causal %in% last.hit)

  preci <- tp / (tp + fp)
  recal <- tp / (tp + fn)

  f1 <- 2 * (preci * recal) / (preci + recal)

  preci <- ifelse(is.na(preci), 0, preci)
  recal <- ifelse(is.na(recal), 0, recal)
  f1 <- ifelse(is.na(f1), 0, f1)

  return(list(f1 = f1, precision = preci, recall = recal))
}

#' F1.LD : Calcul F1 score for a simulation (on effect sizes), take in account the linkage desequilibrum
#'
#' @param beta : effect from a regression
#' @param causal : true causal variable (SNP or CpG)
#' @param nb.hit : number of first hit for calculate the F1 score
#'
#' @return the F1 score, the precision and the recall
#'
#' @export
F1.LD <- function(beta, causal, nb.hit, neighbour.c = 50, neighbour.p = 50) {

  pos <- sort(order(abs(beta), decreasing = T)[1:nb.hit])

  # proche voisin
  pN <- NULL
  for (i in pos) {
    ai <- i
    a1 <- ai - neighbour.c
    a2 <- ai + neighbour.c
    pN <- c(pN, a1:a2)
  }
  pN <- unique(sort(pN))

  ###### Compter le nombre de pique

  NbP <- NULL
  for (j in 1:length(pos)) {
    pos1 <- pos[j]
    pos2 <- pos[j + 1]

    pos50 <- pos1:(pos1 + neighbour.p)

    if ((pos1 %in% pos50) & (pos2 %in% pos50)) {
      NbP <- c(NbP, 0)
    }
    else {NbP <- c(NbP, 1)}
  }

  np <- sum(NbP)

  # TP oui ou non
  tp <- sum(pN %in% causal)
  # FP oui ou non
  fp <- np - tp
  fp <- ifelse(fp < 0, 0, fp)
  # FN oui ou non
  fn <- sum(!(causal %in% pN))

  preci <- tp / (tp + fp)
  recal <- tp / (tp + fn)

  f1 <- 2 * (preci * recal) / (preci + recal)

  preci <- ifelse(is.na(preci), 0, preci)
  recal <- ifelse(is.na(recal), 0, recal)
  f1 <- ifelse(is.na(f1), 0, f1)

  return(list(f1 = f1, precision = preci, recall = recal))
}


#' F1.SCORE.FOR.HIMA : F1 score for HIMA package (for simulation)
#'
#' @param hima.data data of hima
#' @param causal true result of the simulation
#' @param M matrix of methylation
#' @return precision : statistical precision of the test : TP / (TP + FP)
#' @return recall : statistical recall of the test : TP / (TP + FN)
#' @return f1_score : F1 Score of the test : 2(Pre * Rec) / (Pre + Rec)
#' @return length.list : number of pValues lower than "ral" (Bonferroni correction).
#' @return TP : True positive
#' @return FP : False positive
#' @return FN : False negative
#'
#' @details
#'
#' F1 = 2*recall*power/(recall + power)
#' @export
F1.SCORE.FOR.HIMA <- function(hima.data, causal, M = NULL) {
  pos <- as.numeric(gsub("`", "", row.names(him)))
  neg <- (1:ncol(M))[-pos]

  TP <- sum(pos %in% causal)
  FP <- sum(!(pos %in% causal))
  FN <- sum(neg %in% causal)

  pre <- TP / (TP + FP)
  rec <- TP / (TP + FN)

  f1_score = 2*rec*pre/(rec + pre)


  return(list(precision = pre,
              recall = rec,
              f1_score = 2*rec*pre/(rec + pre),
              length.list = length(pos),
              TP = TP,
              FP = FP,
              FN = FN))
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

#' refactor2 : fonction refactor
#'
#' @param data methylation matrix
#' @param k number of K
#' @param covar covariables
#' @param t number of of CpGs for the estimation of latent factor
#' @param numcomp like K
#' @param stdth sites were excluded due to low variance
#'
#'
#' @details
#'
#' see refactor paper
#' @export
refactor2 <- function(data, k, covar = NULL, t = 500, numcomp = NULL, stdth = 0.02) {

  #ranked_filename = paste(out, ".out.rankedlist.txt", sep="")
  #components_filename = paste(out, ".out.components.txt", sep="")

  # print('Starting ReFACTor v1.0...');

  # print('Reading input files...');

  O <- as.matrix(data)
  # sample_id <- O[1, -1] # extract samples ID
  # O <- O[-1,] # remove sample ID from matrix
  # cpgnames <- O[, 1] ## set rownames
  # O <- O[, -1]
  # O = matrix(as.numeric(O),nrow=nrow(O),ncol=ncol(O))
  #
  # print(paste("Excluding sites with low variance (std < ", stdth, ")..."), sep="")
  sds <- apply(t(O), 2, sd)
  m_before <- length(sds)
  include <- which(sds >= stdth)
  O <- O[include,]
  # cpgnames = cpgnames[include]
  # print(paste((m_before - length(which(sds >= stdth))), " sites were excluded due to low variance...", sep=""))

  if (is.null(numcomp) || is.na(numcomp))
  {
    numcomp <- k
  }

  # Adjust the data for the covariates
  if (!is.null(covar))
  {
    covs <- as.matrix(covar)
    # sample_id2 <- covs[, 1]
    # if (!all(sample_id == sample_id2)){
    #   print("ERROR: The order of the samples in the covariates file must be the same as the order in the data file")
    #   quit()
    # }
    # covs <- covs[,-1]
    if (length(covs) > dim(O)[2])
    {
      covs <- matrix(as.numeric(covs), nrow = nrow(covs), ncol = ncol(covs))
    }else{
      covs <- as.numeric(covs)
    }
    O_adj <- O
    for (site in 1:nrow(O))
    {
      model <- lm(O[site,] ~  covs)
      O_adj[site,] <- residuals(model)
    }
    O <- O_adj
  }

  # print('Running a standard PCA...')
  pcs <- prcomp(scale(t(O)))

  coeff <- pcs$rotation
  score <- pcs$x

  # print('Compute a low rank approximation of input data and rank sites...')
  x <- score[,1:k] %*% t(coeff[,1:k])
  An <- scale(t(O), center = T, scale = F)
  Bn <- scale(x, center = T, scale = F)
  An <- t(t(An) * (1/sqrt(apply(An^2, 2, sum))))
  Bn <- t(t(Bn) * (1/sqrt(apply(Bn^2, 2, sum))))


  # Find the distance of each site from its low rank approximation.
  distances <- apply((An-Bn)^2, 2, sum)^0.5
  dsort <- sort(distances, index.return = T)
  ranked_list <- dsort$ix

  # print('Compute ReFACTor components...')
  sites <- ranked_list[1:t]
  pcs <- prcomp(scale(t(O[sites,])))
  first_score <- score[,1:k]
  score <- pcs$x

  #print('Saving a ranked list of the data features...');
  #write(t(cpgnames[ranked_list]),file=ranked_filename,ncol=1)
  #write(t(cbind(ranked_list,cpgnames[ranked_list])),file=ranked_filename,ncol=2)

  #print('Saving the ReFACTor components...');
  #write(t(score[,1:numcomp]), file=components_filename, ncol=numcomp)

  # print('ReFACTor is Done');
  result <- list(refactor_components = score[, 1:numcomp],
                 ranked_list = ranked_list,
                 standard_pca = first_score)
  return(result)

}


#' F1.SCORE.FOR.SPARSE : F1 score for sparse LFMM package (for simulation)
#'
#' @param a alpha effect in mediation (sparse_lfmm$B[,1])
#' @param b beta effect in mediation (sparse_lfmm$B[,2])
#' @param causal true result of the simulation
#' @return precision : statistical precision of the test : TP / (TP + FP)
#' @return recall : statistical recall of the test : TP / (TP + FN)
#' @return f1_score : F1 Score of the test : 2(Pre * Rec) / (Pre + Rec)
#' @return length.list : number of pValues lower than "ral" (Bonferroni correction).
#' @return TP : True positive
#' @return FP : False positive
#' @return FN : False negative
#'
#' @details
#' Select the markers with a * b != 0
#' F1 = 2*recall*power/(recall + power)
#' @export
F1.SCORE.FOR.SPARSE <- function(a, b, causal) {
  pos <- which((a * b) != 0)
  neg <- (1:length(a))[-pos]

  TP <- sum(pos %in% causal)
  FP <- sum(!(pos %in% causal))
  FN <- sum(neg %in% causal)

  pre <- TP / (TP + FP)
  rec <- TP / (TP + FN)

  f1_score = 2*rec*pre/(rec + pre)


  return(list(precision = pre,
              recall = rec,
              f1_score = 2*rec*pre/(rec + pre),
              length.list = length(pos),
              TP = TP,
              FP = FP,
              FN = FN))
}


#' F1.SCORE.FOR.BMA : F1 score for bigmma package (for simulation)
#'
#' @param bma result of the data.org.big function
#' @param causal true result of the simulation
#' @param M matrix of methylation
#' @return precision : statistical precision of the test : TP / (TP + FP)
#' @return recall : statistical recall of the test : TP / (TP + FN)
#' @return f1_score : F1 Score of the test : 2(Pre * Rec) / (Pre + Rec)
#' @return length.list : number of pValues lower than "ral" (Bonferroni correction).
#' @return TP : True positive
#' @return FP : False positive
#' @return FN : False negative
#'
#' @details
#'
#' F1 = 2*recall*power/(recall + power)
#' @export
F1.SCORE.FOR.BMA <- function(bma, causal, M = NULL) {
  pos <- as.numeric(gsub("X", "", summary(bma, only=TRUE)[["mediator"]]))
  neg <- (1:ncol(M))[-pos]

  TP <- sum(pos %in% causal)
  FP <- sum(!(pos %in% causal))
  FN <- sum(neg %in% causal)

  pre <- TP / (TP + FP)
  rec <- TP / (TP + FN)

  f1_score = 2*rec*pre/(rec + pre)


  return(list(precision = pre,
              recall = rec,
              f1_score = 2*rec*pre/(rec + pre),
              length.list = length(pos),
              TP = TP,
              FP = FP,
              FN = FN))
}
