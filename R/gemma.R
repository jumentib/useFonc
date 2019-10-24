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
